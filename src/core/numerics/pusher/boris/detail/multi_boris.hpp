#ifndef PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP

#include "core/utilities/thread_pool.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

enum class MultiBorisMode : std::uint16_t { REF = 0, COPY };

struct MultiBorisOptions
{
    bool use_main_thread = false; // true for perf
};

// ── MultiBoris state struct for AoSPCTS / ModelAccessor ───────────────────────────────

template<typename ModelAccessor, auto _opts = MultiBorisOptions{}>
struct MultiBoris
{
    static constexpr auto opts = _opts;
    static constexpr auto dim  = ModelAccessor::GridLayout_t::dimension;
    using Model_t              = ModelAccessor::Model_t;
    using GridLayout_t         = ModelAccessor::GridLayout_t;
    using ParticleArray_t      = Model_t::particle_array_type;
    using Electromag_t         = Model_t::electromag_type;
    using Vecfield_t           = Electromag_t::vecfield_type;
    using Field_t              = Vecfield_t::field_type;
    using ParticleArray_v      = ParticleArray_t::view_t;

    MultiBoris(double const dt_, ModelAccessor& _accessor, std::function<void(int)> fn_ = {})
        : dt{dt_}
        , accessor{_accessor}
        , fn{fn_}
    {
    }

    double const dt;
    ModelAccessor& accessor;
    std::function<void(int)> fn;

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](double const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};


// ── Per-tile functors (shared between CPU and GPU specialisations) ─────────────────────

template<auto particle_type, auto boris_mode, typename MultiBorisPusherImpl_t>
struct MultiBorisFunctors
{
    static_assert(all_are<ParticleType>(particle_type));

    using GridLayout_t    = MultiBorisPusherImpl_t::GridLayout_t;
    using Particles_t     = MultiBorisPusherImpl_t::Particles_t;
    using Electromag_t    = MultiBorisPusherImpl_t::Electromag_t;
    using Interpolator_t  = MultiBorisPusherImpl_t::Interpolator_t;
    using ParticleArray_v = Particles_t::view_t;
    using Vecfield_t      = Electromag_t::vecfield_type;
    using Field_t         = Vecfield_t::field_type;
    using Tile_vt         = Field_t::value_type::value_type;
    using VecField_vt     = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt   = basic::Electromag<VecField_vt>;

    static constexpr auto dim = GridLayout_t::dimension;
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    MultiBorisFunctors(auto& in, auto& view, auto& pop, auto& parts, auto& em)
        : pps{*parts}
        , electromag{em}
        , dto2m{0.5 * in.dt / pop.mass()}
        , halfdt{in.mesh(view.layout.meshSize(), in.dt)}
    {
    }

    void operator()(auto& in, [[maybe_unused]] auto const i) { on_cpu_tiles(in); }

    void on_cpu_tiles(auto& /*in*/)
    {
        for (std::size_t tileidx = 0; tileidx < pps().size(); ++tileidx)
            one_tile(tileidx);
    }

    // one tile's worth of work — the unit of dispatch for pool.detach_task() in the
    // thread-pooled path (see MultiBorisPusherImpl::move_rest_pooled).
    void one_tile(std::size_t const tile_idx) { per_tile(tile_idx); }

    void per_tile(std::size_t const tile_idx)
    {
        auto& tile         = pps()[tile_idx];
        auto const& layout = electromag.E[0][tile_idx].layout();
        auto const& em     = em_tile(tile_idx);

        auto& parts          = tile();
        auto const tile_cell = pps.local_cell(tile.lower);

        auto constexpr static tracker = [](auto&&... args) {
            return make_particle_tracker<Particles_t::layout_mode, particle_type, dim>(args...);
        };

        for (auto const& bix : parts.local_box())
        {
            auto& cell_particles = parts.particles_(bix);
            for (std::size_t pid = 0; pid < cell_particles.size(); ++pid)
            {
                auto const& old_cell = cell_particles[pid].iCell();
                // tile_cell names the tile PHYSICALLY holding this particle;
                // pps.local_tile_cell would instead give the clamp-owner, which
                // diverges for duplicated level ghost cells
                auto const pt = tracker(old_cell, tile_cell);
                per_particle(cell_particles[pid], layout, pt, pid, em);
            }
        }
    }

    void per_particle_still_in_ghost_box(auto&&... args)
    {
        static constexpr auto alloc_mode             = Particles_t::alloc_mode;
        auto const& [particle, layout, pt, pidx, em] = std::forward_as_tuple(args...);

        {
            Interpolator_t interp;
            boris::accelerate(particle, interp.m2p(particle, em, layout), dto2m);
        }
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

        if constexpr (boris_mode == MultiBorisMode::REF)
            // per-cell buckets are tracked even in the tile ghost layer, so any particle
            // whose cell changed needs registering — domain and level ghost alike. pt was
            // built in per_tile, before advance() ran.
            pps.template move_check<particle_type>(pt, pidx, particle);
    }

    void per_particle(auto&&... args)
    {
        static constexpr auto alloc_mode             = Particles_t::alloc_mode;
        auto const& [particle, layout, pt, pidx, em] = std::forward_as_tuple(args...);

        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

        if constexpr (particle_type == ParticleType::Domain)
            per_particle_still_in_ghost_box(args...);
        else if constexpr (particle_type == ParticleType::LevelGhost)
        {
            if (isIn(particle, pps.ghost_box()))
                per_particle_still_in_ghost_box(args...);
            else if constexpr (boris_mode == MultiBorisMode::REF)
                // left the ghost box on the first half-step — register the departure or
                // the bucket holding it goes stale; move_check resolves it as a deletion
                pps.template move_check<particle_type>(pt, pidx, particle);
        }
        else
        {
            PHARE_ASSERT(false);
        }
    }

    // one tile's worth of the COPY-mode deposit pass. Domain and LevelGhost write into
    // the SAME tile's rhoP/rhoC/flux, so callers must run this tile's Domain and
    // LevelGhost passes in one task (never two concurrent tasks for the same tile_idx)
    // — see MultiBorisPusherImpl::move_cpu_copy_pooled. Different tiles are independent.
    void one_copy_tile(std::size_t const tile_idx, auto& boxings, auto& view, auto& pop)
    {
        auto const& patch_id      = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = pps.box();

        auto& pctile         = pps()[tile_idx];
        auto& cell_particles = pctile();
        bool const is_border = patch_box * grow(pctile, 1) != grow(pctile, 1);
        auto const& layout   = electromag.E[0][tile_idx].layout();
        auto const tile_em   = em_tile(tile_idx);
        auto& rhoP           = pop.particleDensity()[tile_idx];
        auto& rhoC           = pop.chargeDensity()[tile_idx];
        auto F = pop.flux().template as<VecField_vt>([&](auto& c) { return c()[tile_idx](); });
        Interpolator_t interp;

        for (auto& cell : cell_particles()) // all cells, for levelghosts
        {
            for (auto p : cell)
            {
                if constexpr (particle_type == ParticleType::LevelGhost)
                    if (!isIn(p, pps.ghost_box()))
                        continue;

                p.iCell() = boris::advance<AllocatorMode::CPU>(p, halfdt);
                boris::accelerate(p, interp.m2p(p, tile_em, layout), dto2m);
                p.iCell() = boris::advance<AllocatorMode::CPU>(p, halfdt);

                if constexpr (particle_type == ParticleType::Domain)
                {
                    if (!is_border || isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                        interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                }
                else if constexpr (particle_type == ParticleType::LevelGhost)
                {
                    // only this tile's own duplicate (checking the whole patch would
                    // double-count level ghosts duplicated across tile borders)
                    if (isIn(p.iCell(), pctile))
                        interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                }
            }
        }
    }

    void per_copy_of_cpu_tile(auto& boxings, auto& view, auto& pop)
    {
        for (std::size_t tile_idx = 0; tile_idx < pps().size(); ++tile_idx)
            one_copy_tile(tile_idx, boxings, view, pop);
    }

    auto em_tile(auto const tidx)
    {
        return electromag.template as<Electromag_vt>([&](auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx](); });
        });
    }

    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const dto2m;
    std::array<double, dim> halfdt;
};


// AoSPCTS only: tiles of per-cell AoS particles, tiled fields. Particle iteration relies
// on each tile containing a PerCellParticles container (flat iterators over all particles
// in that tile), so this isn't generic over LayoutMode - `layout` is only kept as a
// template parameter so callers can name it the same way as other per-layout impls.
template<LayoutMode layout, AllocatorMode alloc, typename GridLayout, typename Particles,
         typename Electromag, typename Interpolator>
struct MultiBorisPusherImpl
{
    static_assert(layout == LayoutMode::AoSPCTS, "MultiBorisPusherImpl only supports AoSPCTS");

    using This
        = MultiBorisPusherImpl<layout, alloc, GridLayout, Particles, Electromag, Interpolator>;

    static constexpr auto dim = GridLayout::dimension;

    using GridLayout_t   = GridLayout;
    using Particles_t    = Particles;
    using Electromag_t   = Electromag;
    using Interpolator_t = Interpolator;

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;

    template<auto type>
    static void sync_particles(auto& particles)
    {
        particles.template on_moved<type>();
    }

    template<auto mode, typename ModelAccessor>
    static void move_rest(MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            Functors<ParticleType::Domain, mode>{in, view, pop, domain, em}(in, i);

            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();
            Functors<ParticleType::LevelGhost, mode>{in, view, pop, level_ghost, em}(in, i);
        }
    }

    template<auto mode, typename ModelAccessor>
    static void move_cpu_copy(MultiBoris<ModelAccessor>& in, auto& boxings, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
            parts.reset_views();
            Functors<particle_type, mode> fns{in, view, pop, parts, em};
            fns.per_copy_of_cpu_tile(boxings, view, pop);
        };

        for (auto& pop : ions)
        {
            per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
            per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
        }
    }

    // ── Thread-pooled variants: level 1 (patch -> pool) is dispatched by the caller
    // (move_pooled), which hands us the specific pool already assigned to patch `i`.
    // Here we do level 2: one task per tile, submitted to that SAME pool so all of its
    // threads_per_pool threads work this patch's tiles at once. These must be called
    // from the dispatching thread itself, never from inside a task already running on
    // `pool` — detach_task-then-wait on one's own pool deadlocks.

    template<auto mode, typename ModelAccessor>
    static void move_rest_pooled(auto& pool, MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();

            using DomainFn = Functors<ParticleType::Domain, mode>;
            using GhostFn  = Functors<ParticleType::LevelGhost, mode>;

            // shared_ptr keeps the (self-contained: own pps view + electromag copy)
            // Functors alive for every deferred per-tile task below.
            auto domain_fn = std::make_shared<DomainFn>(in, view, pop, domain, em);
            auto ghost_fn  = std::make_shared<GhostFn>(in, view, pop, level_ghost, em);

            auto const n_tiles = domain_fn->pps().size();
            assert(n_tiles == ghost_fn->pps().size());

            for (std::size_t t = 0; t < n_tiles; ++t)
                pool.detach_task([domain_fn, ghost_fn, t] {
                    domain_fn->one_tile(t);
                    ghost_fn->one_tile(t);
                });
        }
    }

    template<auto mode, typename ModelAccessor>
    static void move_cpu_copy_pooled(auto& pool, MultiBoris<ModelAccessor>& in, auto& boxings,
                                     auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();

            using DomainFn = Functors<ParticleType::Domain, mode>;
            using GhostFn  = Functors<ParticleType::LevelGhost, mode>;

            auto domain_fn = std::make_shared<DomainFn>(in, view, pop, domain, em);
            auto ghost_fn  = std::make_shared<GhostFn>(in, view, pop, level_ghost, em);

            auto const n_tiles = domain_fn->pps().size();
            assert(n_tiles == ghost_fn->pps().size());

            // view/pop are cheap handles, safe to capture by copy. Domain and LevelGhost
            // for the SAME tile must run in the SAME task (shared rhoP/rhoC/flux);
            // different tiles are independent.
            for (std::size_t t = 0; t < n_tiles; ++t)
                pool.detach_task([domain_fn, ghost_fn, t, view, pop, &boxings]() mutable {
                    domain_fn->one_copy_tile(t, boxings, view, pop);
                    ghost_fn->one_copy_tile(t, boxings, view, pop);
                });
        }
    }

    template<typename ModelAccessor>
    static void sync_ref(MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view      = in.accessor[i];
        auto [ions, _] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            sync_particles<ParticleType::Domain>(domain);

            auto& level_ghost = pop.levelGhostParticles();
            sync_particles<ParticleType::LevelGhost>(level_ghost);
        }
    }

    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelAccessor>
    static void move(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        if constexpr (MultiBorisOptions{}.use_main_thread)
            move_sequential<mode>(in, boxings);
        else
            move_pooled<mode>(in, boxings);
    }

    // No thread pools at all: patches then tiles, both fully sequential on the caller's
    // thread. Used when MultiBorisOptions::use_main_thread is set (perf comparisons /
    // environments where spinning up pools isn't worth it).
    template<MultiBorisMode mode, typename ModelAccessor>
    static void move_sequential(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        static constexpr auto copy   = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu = Particles_t::alloc_mode == AllocatorMode::CPU;

        for (std::size_t i = 0; i < in.accessor.size(); ++i)
        {
            if constexpr (copy and is_cpu)
                move_cpu_copy<mode>(in, boxings, i);
            else
                move_rest<mode>(in, i);

            if constexpr (not copy)
                sync_ref(in, i);
        }
    }

    // General case: two levels of parallelism. Each patch is assigned one pool
    // (round-robin over whichever is next free), and within that pool every tile of
    // the patch is its own task -- so with N pools of M threads, up to N patches and
    // N*M tiles are in flight simultaneously. Dispatch happens entirely from this
    // (calling) thread so no pool ever detach_task-then-waits on itself.
    template<MultiBorisMode mode, typename ModelAccessor>
    static void move_pooled(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        static constexpr auto copy   = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu = Particles_t::alloc_mode == AllocatorMode::CPU;

        auto& TP = ThreadPool::INSTANCE();

        for (std::size_t i = 0; i < in.accessor.size(); ++i)
        {
            auto& pool = TP.get_pool(TP.first_ready_idx());
            if constexpr (copy and is_cpu)
                move_cpu_copy_pooled<mode>(pool, in, boxings, i);
            else
                move_rest_pooled<mode>(pool, in, i);
        }
        TP.sync(); // every patch's tiles, across every pool, are done

        if constexpr (not copy)
        {
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
            {
                auto& pool = TP.get_pool(TP.first_ready_idx());
                pool.detach_task([&in, i] { sync_ref(in, i); });
            }
            TP.sync();
        }
    }
};

} // namespace PHARE::core


#endif
