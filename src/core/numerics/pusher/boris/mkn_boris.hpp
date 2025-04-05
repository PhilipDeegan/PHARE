#ifndef PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP


#if PHARE_HAVE_MKN_GPU

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
auto static const multi_boris_threads = get_env_as("PHARE_ASYNC_THREADS", std::size_t{5});

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
    using ParticleArray_v     = typename ParticleArray_t::view_t;
    using Box_t               = Box<std::uint32_t, dim>;
    using Boxes_t             = std::vector<Box_t>;
    using Particles_ptrs      = std::vector<ParticleArray_t*>;
    using StreamLauncher      = gpu::ThreadedBoxStreamLauncher<Boxes_t, Particles_ptrs>;


    static auto _particles(ModelViews& views)
    {
        std::vector<ParticleArray_t*> ptrs;
        auto add = [&](auto& ps) {
            // if (ps.size())
            ptrs.emplace_back(&ps);
        };
        auto all = [&](auto&... ps) { (add(ps), ...); };

        for (auto& view : views)
            for (auto& pop : view.ions)
                all(pop.domainParticles(), pop.patchGhostParticles(), pop.levelGhostParticles());

        return ptrs;
    }

    MultiBoris(double const dt_, ModelViews& _views)
        : dt{dt_}
        , views{_views}
        , particles{_particles(views)}
    {
        for (auto& view : views)
        {
            auto each = [&](auto const& type, auto& ps, auto& pop) mutable {
                // if (!ps.size())
                //     return;
                particle_type.emplace_back(type);
                auto& domainview = pviews.emplace_back(*ps);
                boxes.emplace_back(domainview.local_box());
                layouts.emplace_back(view.layout);
                ems.emplace_back(view.electromag);
                dto2ms.emplace_back(0.5 * dt / pop.mass());
                halfdt.emplace_back(mesh(view.layout.meshSize(), dt));
                rhos.emplace_back(pop.density());
                fluxes.emplace_back(pop.flux());
            };

            for (auto& pop : *view.ions)
            {
                each(0, pop.domainParticles(), pop);
                each(1, pop.patchGhostParticles(), pop);
                each(2, pop.levelGhostParticles(), pop);
            }
        }
    }

    void reset() {}

    double dt;
    ModelViews& views;

    // on host
    Particles_ptrs particles;
    std::vector<std::uint16_t> particle_type;
    Boxes_t boxes;

    // on device
    gpu::Vec_t<ParticleArray_v> pviews;
    gpu::Vec_t<GridLayout_t> layouts;
    gpu::Vec_t<Electromag_t> ems;
    gpu::Vec_t<Field_t> rhos;
    gpu::Vec_t<Vecfield_t> fluxes;
    gpu::Vec_t<std::array<double, dim>> halfdt;
    gpu::Vec_t<double> dto2ms;

    StreamLauncher streamer{particles, boxes, detail::multi_boris_threads};

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

public:
    template<typename ModelViews>
    static void move(MultiBoris<ModelViews>& in)
    {
        using ParticleArray_v            = typename MultiBoris<ModelViews>::ParticleArray_v;
        static constexpr auto alloc_mode = Particles_t::alloc_mode;

        auto const emps         = in.ems.data();
        auto const dto2mspp     = in.dto2ms.data();
        auto const layoutps     = in.layouts.data();
        auto const halfDtOverDl = in.halfdt.data();
        auto pps                = in.pviews.data();

        auto per_particle = [=] _PHARE_ALL_FN_(auto&&... args) mutable {
            auto const& [particle, view, tile_cell, pidx, i] = std::forward_as_tuple(args...);

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
            {
                Interpolator interp;
                boris::accelerate(particle, interp.m2p(particle, emps[i], layoutps[i]),
                                  dto2mspp[i]);
            }
            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

            auto const new_cell   = view.local_tile_cell(particle.iCell());
            bool const moved_tile = !array_equals(new_cell, tile_cell);
            if (moved_tile)
                view.icell_changer(particle, tile_cell, pidx, particle.iCell());
        };

        auto per_any_particle = [=] _PHARE_DEV_FN_(auto& particles, auto&&... args) mutable {
            auto const& pidx = std::get<2>(std::forward_as_tuple(args...));

            using enum LayoutMode;
            if constexpr (any_in(ParticleArray_v::layout_mode, SoA, SoAPC, SoATS))
                per_particle(detail::SoAZipParticle{particles, pidx}, args...);
            else
                per_particle(particles[pidx], args...);
        };

        auto per_tile = [=] _PHARE_DEV_FN_(auto const& i) mutable {
            auto& view      = pps[i];
            auto& tile      = view()[blockIdx.x];
            auto& parts     = tile();
            auto const tidx = threadIdx.x;
            auto const ws   = kernel::warp_size();
            auto const each = parts.size() / ws;

            std::size_t pid = 0;
            for (; pid < each; ++pid)
                per_any_particle(parts, view, view.local_cell(tile.lower), pid * ws + tidx, i);
            if (tidx < parts.size() - (ws * each))
                per_any_particle(parts, view, view.local_cell(tile.lower), pid * ws + tidx, i);
        };


        auto& streamer = in.streamer;
        auto ip        = &in; // used in lambdas, copy address! NO REF!

        if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            streamer.async_dev_tiles([=] _PHARE_DEV_FN_(auto const i) mutable { per_tile(i); });

            streamer.host([ip = ip](auto const i) mutable {
                if (ip->particles[i]->size() == 0)
                    return;

                constexpr static std::uint32_t PHASE = 1;
                if (ip->particle_type[i] == 0) //  :(
                    This::template sync<PHASE, ParticleType::Domain>(ip, i);
                else if (any_in(ip->particle_type[i], 1, 2))
                    This::template sync<PHASE, ParticleType::Ghost>(ip, i);
                else
                    throw std::runtime_error("No impl");
            });
        }
        else
        {
            // do cpu
            streamer.host([=](auto const i) mutable {
                if (ip->particles[i]->size() == 0)
                    return;

                auto& particles_ts = *ip->particles[i];

                for (auto& tile : particles_ts())
                {
                    auto const locell  = particles_ts.local_cell(tile.lower);
                    auto const& nparts = particles_ts.size(locell);
                    for (std::size_t pidx = 0; pidx < nparts; ++pidx)
                        per_any_particle(tile(), particles_ts, locell, pidx, i);
                }
            });

            // separation for profiling the same as gpu
            streamer.host([=](auto const i) mutable {
                constexpr static std::uint32_t PHASE = 1;

                auto& particles_ts = *ip->particles[i];
                if (ip->particle_type[i] == 0)
                    ParticleArrayService::sync<PHASE, ParticleType::Domain>(particles_ts);
                else if (any_in(ip->particle_type[i], 1, 2))
                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(particles_ts);
                else
                    throw std::runtime_error("No impl");
            });
        }
    }

    template<auto phase, auto type>
    static void sync(auto p, auto const& i)
    {
        assert(i < p->pviews.size() and i < p->streamer.streams.size());
        p->pviews[i].template sync<phase, type>(p->streamer.streams[i]);
        ParticleArrayService::sync<phase, type>(*p->particles[i]);
    }
};


} // namespace PHARE::core

#endif // PHARE_HAVE_MKN_GPU

#endif /*PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP*/
