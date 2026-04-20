// IWYU pragma: private, include "core/numerics/pusher/multi_boris.hpp"

#ifndef PHARE_CORE_PUSHER_BORIS_MKN_BORIS_CAPTURE_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_BORIS_CAPTURE_HPP


#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"



#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"

// #include "core/numerics/pusher/boris/move.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_service.hpp"


#include <cmath>
#include <cassert>
#include <cstddef>


namespace PHARE::core::detail
{
auto static const multi_boris_threads = get_env_as("PHARE_ASYNC_THREADS", std::size_t{1});

} // namespace PHARE::core::detail

namespace PHARE::core
{

// !!naive!!
template<typename ModelViews>
struct MultiBoris
{
    using ModelView    = ModelViews::value_type;
    using GridLayout_t = ModelView::GridLayout_t;

    using ParticleArray_t     = ModelView::ParticleArray_t;
    static constexpr auto dim = ParticleArray_t::dimension;
    using Electromag_t        = ModelView::Electromag_t;
    using Vecfield_t          = Electromag_t::vecfield_type;
    using Field_t             = Vecfield_t::field_type;
    using ParticleArray_v     = ParticleArray_t::view_t;
    using Box_t               = Box<std::uint32_t, dim>;
    using Boxes_t             = std::vector<Box_t>;
    using Particles_ptrs      = std::vector<ParticleArray_t*>;
    using StreamLauncher      = gpu::ThreadedStreamLauncher<ModelViews>;


    MultiBoris(double const dt_, ModelViews& _views)
        : dt{dt_}
        , views{_views}
    {
    }

    void reset() {}

    double const dt;
    ModelViews& views;
    StreamLauncher streamer{views, detail::multi_boris_threads};

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};




template<typename GridLayout, typename Particles, typename Electromag, typename Interpolator>
class MultiBorisPusher
{
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;
    using This                = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator>;


    // using Electromag_t        = ModelView::Electromag_t;
    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Electromag_vt = basic::Electromag<basic::TensorField<typename Field_t::Super>>;

public:
    template<typename ModelViews>
    static void move(MultiBoris<ModelViews>& in)
    {
        static constexpr auto alloc_mode = Particles_t::alloc_mode;

        auto& streamer = in.streamer;
        auto ip        = &in; // used in lambdas, copy address! NO REF!


        for (std::size_t i = 0; i < in.streamer.datas.size(); ++i)
        {
            auto& view         = in.views[i];
            auto const& layout = view.layout;
            auto const halfdt  = in.mesh(layout.meshSize(), in.dt);
            auto const em      = view.electromag.template as<Electromag_vt>([](auto const& vf) {
                return core::for_N<3, core::for_N_R_mode::make_array>(
                    [&](auto i) { return *vf[i]; });
            });

            auto per_pop = [&](auto& pop) {
                auto const dto2m = 0.5 * in.dt / pop.mass();

                auto per_particles = [&](auto& particle_array) {
                    if (particle_array.size() == 0)
                        return;

                    auto pps                     = *particle_array;
                    [[maybe_unused]] auto ps_ptr = &particle_array; // only for CPU

                    auto per_tile = [=] _PHARE_ALL_FN_(auto const& tile_picker) mutable {
                        auto per_particle = [=] _PHARE_ALL_FN_(auto&&... args) mutable {
                            auto const& [particle, tile_cell, pidx]
                                = std::forward_as_tuple(args...);

                            particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

                            {
                                Interpolator interp;
                                boris::accelerate(particle, interp.m2p(particle, em, layout),
                                                  dto2m);
                            }
                            particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

                            auto const new_cell   = pps.local_tile_cell(particle.iCell());
                            bool const moved_tile = !array_equals(new_cell, tile_cell);
                            if (moved_tile)
                            {
                                if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
                                    pps.icell_changer(particle, tile_cell, pidx, particle.iCell());
                                else // otherwise just push back
                                    ps_ptr->icell_changer(particle, tile_cell, pidx,
                                                          particle.iCell());
                            }
                        };

                        auto per_any_particle = [=] _PHARE_ALL_FN_(auto& particles,
                                                                   auto&&... args) mutable {
                            auto const& pidx = std::get<1>(std::forward_as_tuple(args...));
#if PHARE_HAVE_THRUST
                            using enum LayoutMode;
                            if constexpr (any_in(Particles_t::layout_mode, SoA, SoAPC, SoATS))
                                per_particle(detail::SoAZipParticle{particles, pidx}, args...);
                            else
#endif
                                per_particle(particles[pidx], args...);
                        };

                        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
                        auto& tile                           = *tileptr;
                        auto& parts                          = tile();

                        auto const each      = pps()[tile_idx]().size() / ws; // SPAN SIZE!
                        auto const tile_cell = pps.local_cell(tile.lower);

                        std::size_t pid = 0;
                        for (; pid < each; ++pid)
                            per_any_particle(parts, tile_cell, pid * ws + tidx);
                        if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
                            if (tidx < parts.size() - (ws * each))
                                per_any_particle(parts, tile_cell, pid * ws + tidx);
                    };

                    if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
                    {
#if !PHARE_HAVE_MKN_GPU
                        static_assert(false && "NEEDS GPU IMPL!");
#else
                        using Launcher = gpu::ChunkLauncher<false>;
                        Launcher launcher{1, 0};
                        launcher.b.x = kernel::warp_size();
                        launcher.g.x = particle_array().size();

                        auto const tile_picker = [=] __device__() {
                            // SPAN TILE!
                            return std::make_tuple(1, &pps()[blockIdx.x], threadIdx.x,
                                                   kernel::warp_size());
                        };

                        launcher.stream(in.streamer.streams[i],
                                        [=] __device__() mutable { per_tile(tile_picker); });
#endif // PHARE_HAVE_MKN_GPU
                    }
                    else
                    {
                        std::size_t tileidx    = 0;
                        auto const tile_picker = [&]() { // REF!
                            // NOT SPAN TILE!
                            return std::make_tuple(tileidx, &(*ps_ptr)()[tileidx], 0, 1);
                        };

                        for (; tileidx < particle_array().size(); ++tileidx)
                            per_tile(tile_picker);
                    }
                };

                per_particles(pop.domainParticles());
                // per_particles(pop.patchGhostParticles());
                per_particles(pop.levelGhostParticles());
            };

            for (auto& pop : view.ions)
                per_pop(pop);
        }



        if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            streamer.host([&](auto const i) mutable {
                constexpr static std::uint32_t PHASE = 1;

                auto& view = in.views[i];

                for (auto& pop : view.ions)
                {
                    ParticleArrayService::sync<PHASE, ParticleType::Domain>(pop.domainParticles());
                    (*pop.domainParticles())
                        .template sync<PHASE, ParticleType::Domain>(in.streamer.streams[i]);

                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(
                        pop.patchGhostParticles());
                    (*pop.patchGhostParticles())
                        .template sync<PHASE, ParticleType::Ghost>(in.streamer.streams[i]);

                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(
                        pop.levelGhostParticles());
                    (*pop.levelGhostParticles())
                        .template sync<PHASE, ParticleType::Ghost>(in.streamer.streams[i]);
                }
            });
        }
        else
        {
            streamer.host([&](auto const i) mutable {
                constexpr static std::uint32_t PHASE = 1;

                auto& view = in.views[i];

                for (auto& pop : view.ions)
                {
                    ParticleArrayService::sync<PHASE, ParticleType::Domain>(pop.domainParticles());
                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(
                        pop.patchGhostParticles());
                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(
                        pop.levelGhostParticles());
                };
            });
        }
    }
};

} // namespace PHARE::core


#endif /*PHARE_CORE_PUSHER_BORIS_MKN_BORIS_CAPTURE_HPP*/
