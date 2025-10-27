// IWYU pragma: private, include "core/numerics/pusher/multi_boris.hpp"

#ifndef PHARE_CORE_PUSHER_BORIS_TILE_BORIS_CAPTURE_HPP
#define PHARE_CORE_PUSHER_BORIS_TILE_BORIS_CAPTURE_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"

#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"

// #include "core/numerics/pusher/boris/move.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_service.hpp"

#include "core/utilities/thread_pool.hpp" // defaults to 1 thread, setup during static init!

#include <cmath>
#include <memory>
#include <cassert>
#include <cstddef>


namespace PHARE::core::detail
{
auto static const multi_boris_threads = get_env_as("PHARE_ASYNC_THREADS", std::size_t{1});

} // namespace PHARE::core::detail

namespace PHARE::core
{

enum class MultiBorisMode : std::uint16_t { REF = 0, COPY };

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
    using PerTileParticles_t  = ParticleArray_t::per_tile_particles;
    using Vec_t
        = MinimizingVector<std::shared_ptr<ParticleArray_t>, ParticleArray_t::alloc_mode, 1>;


    MultiBoris(double const dt_, ModelViews& _views, std::function<void(int)> fn_ = {})
        : dt{dt_}
        , views{_views}
        , fn{fn_}
    {
    }

    void reset() {}

    double const dt;
    ModelViews& views;
    std::function<void(int)> fn;
    StreamLauncher streamer{views, detail::multi_boris_threads};

    static inline Vec_t domains; // back up for first push
    static inline Vec_t levelGhosts;

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](double const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};


template<auto particle_type, auto boris_mode, typename MultiBorisPusher_t>
struct MultiBorisFunctors
{
    static_assert(all_are<ParticleType>(particle_type));

    using GridLayout_t    = MultiBorisPusher_t::GridLayout_t;
    using Particles_t     = MultiBorisPusher_t::Particles_t;
    using Electromag_t    = MultiBorisPusher_t::Electromag_t;
    using Interpolator_t  = MultiBorisPusher_t::Interpolator_t;
    using ParticleArray_v = Particles_t::view_t;

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    static constexpr auto dim = GridLayout_t::dimension;
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    MultiBorisFunctors(auto& in, auto& view, auto& pop, auto& parts)
        : pps{*parts}
        , electromag{*view.electromag}
        , dto2m{0.5 * in.dt / pop.mass()}
        , halfdt{in.mesh(view.layout.meshSize(), in.dt)}
    {
        check_particles(parts);
        check_particles_views(parts);
    }

    void check(auto const& particle) _PHARE_ALL_FN_
    {
        PHARE_DEBUG_DO({
            for (std::size_t i = 0; i < GridLayout_t::dimension; ++i)
            {
                assert(std::abs(particle.iCell()[i]) < 1000);
                assert(not std::isnan(particle.delta()[i]));
            }
        })
    }

    void operator()(auto& in, auto const i)
    {
        if constexpr (Particles_t::alloc_mode == AllocatorMode::CPU)
            on_cpu_tiles(in);
        else if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_gpu_tiles(in);
        else
            throw std::runtime_error("Unknown mode");
    }


    void on_gpu_tiles(auto& in)
    {
#if !PHARE_HAVE_GPU
        throw std::runtime_error("NEEDS GPU IMPL!");
#else
        using Launcher = gpu::ChunkLauncher<false>;
        Launcher launcher{1, 0};
        launcher.b.x           = kernel::warp_size();
        launcher.g.x           = pps().size();
        auto const tile_picker = [pps = pps] __device__() {
            return std::make_tuple(blockIdx.x, &pps()[blockIdx.x], threadIdx.x,
                                   kernel::warp_size());
        };

        launcher.stream(in.streamer.streams[i],
                        [=, self = *this] __device__() mutable { self.per_tile(tile_picker); });
#endif // PHARE_HAVE_MKN_GPU
    }

    void on_cpu_tiles(auto& /*in*/)
    {
        std::size_t tileidx    = 0;
        auto const tile_picker = [&]() { return std::make_tuple(tileidx, &pps()[tileidx], 0, 1); };

        for (; tileidx < pps().size(); ++tileidx)
            per_tile(tile_picker);
    }

    void per_tile(auto const& tile_picker) _PHARE_ALL_FN_
    {
        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
        auto& tile                           = *tileptr;
        auto const& layout                   = electromag.E[0][tile_idx].layout(); // any per tile
        auto const& em                       = em_tile(tile_idx);

        using enum LayoutMode;
        if constexpr (boris_mode == MultiBorisMode::REF and any_in(Particles_t::layout_mode))
        {
            aos_cpu_mapped_ref(tile, layout, em);
        }
        else
        {
            auto& parts          = tile();
            auto const each      = pps()[tile_idx]().size() / ws; // SPAN SIZE!
            auto const tile_cell = pps.local_cell(tile.lower);

            std::size_t pid = 0;
            for (; pid < each; ++pid)
                per_any_particle(parts, layout, tile_cell, pid * ws + tidx, em);
            if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
                if (tidx < parts.size() - (ws * each))
                    per_any_particle(parts, layout, tile_cell, pid * ws + tidx, em);
        }
    }
    void per_any_particle(auto& particles, auto&&... args) _PHARE_ALL_FN_
    {
        auto const& pidx = std::get<2>(std::forward_as_tuple(args...));
#if PHARE_HAVE_THRUST
        using enum LayoutMode;
        if constexpr (any_in(Particles_t::layout_mode, SoA, SoAPC, SoATS))
            per_particle(SoAZipParticle{particles, pidx}, args...);
        else
#endif
            per_particle(particles[pidx], args...);
    }


    void per_particle_still_in_ghost_box(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode                    = Particles_t::alloc_mode;
        auto const& [particle, layout, tile_cell, pidx, em] = std::forward_as_tuple(args...);

        {
            Interpolator_t interp;
            boris::accelerate(particle, interp.m2p(particle, em, layout), dto2m);
        }
        check(particle);
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);
        check(particle);

        if constexpr (boris_mode == MultiBorisMode::REF and particle_type == ParticleType::Domain)
            if (isIn(particle, pps.box()))
            {
                auto const new_cell   = pps.local_tile_cell(particle.iCell());
                bool const moved_tile = !array_equals(new_cell, tile_cell);
                if (moved_tile)
                    pps.icell_changer(tile_cell, pidx, particle.iCell());
            }
        check(particle);
    }


    void per_particle(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode                    = Particles_t::alloc_mode;
        auto const& [particle, layout, tile_cell, pidx, em] = std::forward_as_tuple(args...);

        check(particle);
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);
        check(particle);

        if constexpr (particle_type == ParticleType::Domain)
            per_particle_still_in_ghost_box(args...);
        else if constexpr (particle_type == ParticleType::LevelGhost)
        {
            if (isIn(particle, pps.ghost_box()))
                per_particle_still_in_ghost_box(args...);
        }
        else
        {
            PHARE_ASSERT(false);
        }
    }


    void aos_cpu_mapped_ref_interp(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode                         = Particles_t::alloc_mode;
        auto const& [og_icell, particle, cidx, tile, layout, em] = std::forward_as_tuple(args...);

        {
            Interpolator_t interp;
            boris::accelerate(particle, interp.m2p(particle, em, layout), dto2m);
        }
        check(particle);
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);
        check(particle);

        if constexpr (boris_mode == MultiBorisMode::REF)
        {
            auto const old_cell   = pps.local_cell(og_icell);
            auto const new_cell   = pps.local_cell(particle.iCell());
            bool const moved_cell = !array_equals(new_cell, old_cell);
            if (moved_cell)
                pps.template icell_changer<particle_type>(old_cell, cidx, particle.iCell());
        }
    }
    void aos_cpu_mapped_ref(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode = Particles_t::alloc_mode;
        auto const& [tile, layout, em]   = std::forward_as_tuple(args...);
        auto const ghost_box             = grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts());

        for (auto const& bix : ghost_box)
        {
            auto const& cell_dexes = tile().map(bix);
            for (std::size_t cidx = 0; cidx < cell_dexes.size(); ++cidx)
            {
                auto const& pidx    = cell_dexes[cidx];
                auto& particle      = tile()[pidx];
                auto const og_icell = particle.iCell();
                particle.iCell()    = boris::advance<alloc_mode>(particle, halfdt);

                if constexpr (particle_type == ParticleType::Domain)
                    aos_cpu_mapped_ref_interp(og_icell, particle, cidx, tile, layout, em);
                else if constexpr (particle_type == ParticleType::LevelGhost)
                {
                    if (isIn(particle, pps.ghost_box()))
                        aos_cpu_mapped_ref_interp(og_icell, particle, cidx, tile, layout, em);
                }
                else
                {
                    PHARE_ASSERT(false);
                }
            }
        }
    }


    auto em_tile(auto const tidx)
    {
        return electromag.template as<Electromag_vt>([&] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx]; });
        });
    }

    void per_copy_of_cpu_tile(auto& boxings, auto& view, auto& pop, auto& parts)
    {
        auto const& patch_id = view.patchID();
        assert(boxings.count(patch_id));
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = parts.box();
        auto const is_border_     = [&](auto const& tile) {
            auto const tile_grow_box = grow(tile, 1);
            return patch_box * tile_grow_box != tile_grow_box;
        };

        for (std::size_t tile_idx = 0; tile_idx < parts().size(); ++tile_idx)
        {
            auto& tile           = parts()[tile_idx];
            auto const tile_cell = parts.local_cell(tile.lower);
            auto const& layout   = electromag.E[0][tile_idx].layout(); // any per tile
            auto const em        = em_tile(tile_idx);
            auto& rhoP           = pop.particleDensity()[tile_idx];
            auto& rhoC           = pop.chargeDensity()[tile_idx];
            auto F = pop.flux().template as<VecField_vt>([&](auto& c) { return c()[tile_idx]; });
            Interpolator_t interp;

            auto on_tile = [&]<bool border>() {
                for (auto p : tile()) // copy particle!
                {
                    per_particle(p, layout, tile_cell, 0, em);
                    if constexpr (!border)
                        interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    else
                    {
                        if (isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                            interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    }
                }
            };

            if (is_border_(tile))
                on_tile.template operator()<true>();
            else
                on_tile.template operator()<false>();
        }
    }


    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const dto2m;
    std::array<double, dim> halfdt;
};

template<typename GridLayout, typename Particles, typename Electromag, typename Interpolator>
class MultiBorisPusher
{
public:
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;
    using This                = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator>;

    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;


    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelViews>
    static void move(MultiBoris<ModelViews>& in, auto const& boxings)
    {
        static constexpr auto alloc_mode = Particles_t::alloc_mode;
        static constexpr auto copy       = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu     = Particles_t::alloc_mode == AllocatorMode::CPU;

        if constexpr (copy and !is_cpu)
        {
            MultiBoris<ModelViews>::domains.get(in.views.size());
            MultiBoris<ModelViews>::levelGhosts.get(in.views.size());
        }

        in.streamer.host([&](auto const i) mutable {
            if constexpr (copy and !is_cpu)
                make_copy(in, i);

            if constexpr (copy and is_cpu)
                move_cpu_copy<mode>(in, boxings, i);

            else
                move_rest<mode>(in, i);
        });

        in.streamer.host([&](auto const i) mutable {
            if constexpr (not copy)
                sync_ref(in, i);
        });
    }


    template<auto mode, typename ModelViews>
    static void move_cpu_copy(MultiBoris<ModelViews>& in, auto& boxings, auto const i)
    {
        auto& view = in.views[i];

        auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
            parts.sync();
            Functors<particle_type, mode> fns{in, view, pop, parts};
            fns.per_copy_of_cpu_tile(boxings, view, pop, parts);
        };

        for (auto& pop : view.ions)
        {
            per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
            per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
        }
    }


    template<typename ModelViews>
    static void make_copy(MultiBoris<ModelViews>& in, auto const i)
    {
        auto& view = in.views[i];

        for (std::size_t j = 0; j < view.ions.size(); ++j)
        {
            auto const popidx = view.ions.size() * i + j;
            auto& pop         = view.ions[j];

            pop.domainParticles().sync();
            MultiBoris<ModelViews>::domains[popidx]
                = std::make_shared<Particles_t>(pop.domainParticles());

            pop.levelGhostParticles().sync();
            MultiBoris<ModelViews>::levelGhosts[popidx]
                = std::make_shared<Particles_t>(pop.levelGhostParticles());
        }
    }


    template<auto mode, typename ModelViews>
    static void move_rest(MultiBoris<ModelViews>& in, auto const i)
    {
        static constexpr auto copy = mode == MultiBorisMode::COPY;

        auto& view = in.views[i];

        for (std::size_t j = 0; j < view.ions.size(); ++j)
        {
            auto const popidx = view.ions.size() * i + j;
            auto& pop         = view.ions[j];

            auto& domain = copy ? *MultiBoris<ModelViews>::domains[popidx] : pop.domainParticles();
            domain.reset_views();
            Functors<ParticleType::Domain, mode>{in, view, pop, domain}(in, i);

            auto& level_ghost
                = copy ? *MultiBoris<ModelViews>::levelGhosts[popidx] : pop.levelGhostParticles();
            level_ghost.reset_views();
            Functors<ParticleType::LevelGhost, mode>{in, view, pop, level_ghost}(in, i);
        }
    }


    template<typename ModelViews>
    static void sync_ref(MultiBoris<ModelViews>& in, auto const i)
    {
        constexpr static std::uint32_t PHASE = 1;

        if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            in.streamer.streams[i].sync();

        auto& view = in.views[i];

        for (std::size_t j = 0; j < view.ions.size(); ++j)
        {
            auto& pop         = view.ions[j];
            auto const popidx = view.ions.size() * i + j;

            auto& domain = pop.domainParticles();
            sync_aos_ts<ParticleType::Domain>(domain, in.streamer.streams[i]);

            auto& level_ghost = pop.levelGhostParticles();
            sync_aos_ts<ParticleType::LevelGhost>(level_ghost, in.streamer.streams[i]);
        }
    }
};

} // namespace PHARE::core


#endif /*PHARE_CORE_PUSHER_BORIS_TILE_BORIS_CAPTURE_HPP*/
