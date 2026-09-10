#ifndef PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP

#include "core/utilities/span.hpp"
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

// Primary (incomplete) template — specialisations provide the implementation per
// (LayoutMode, AllocatorMode) pair.  Add a new #include below for each new layout.
template<LayoutMode layout, AllocatorMode alloc, typename GridLayout, typename Particles,
         typename Electromag, typename Interpolator>
struct MultiBorisPusherImpl;


// ── MultiBoris state struct for AoSTS / ModelAccessor ─────────────────────────────────

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
    // thread-pooled path (see MultiBorisPusherImplBase::move_rest_pooled).
    void one_tile(std::size_t const tile_idx)
    {
        auto const tile_picker
            = [&]() { return std::make_tuple(tile_idx, &pps()[tile_idx], 0, 1); };
        per_tile(tile_picker);
    }

    void per_tile(auto const& tile_picker)
    {
        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
        auto& tile                           = *tileptr;
        auto const& layout                   = electromag.E[0][tile_idx].layout();
        auto const& em                       = em_tile(tile_idx);

        using enum LayoutMode;

        auto& parts          = tile();
        auto const tile_cell = pps.local_cell(tile.lower);

        auto constexpr static tracker = [](auto&&... args) {
            return make_particle_tracker<Particles_t::layout_mode, particle_type, dim>(args...);
        };

        if constexpr (Particles_t::layout_mode == AoSPCTS)
        {
            for (auto const& bix : parts.local_box())
            {
                auto& cell_particles = parts.particles_(bix);
                for (std::size_t pid = 0; pid < cell_particles.size(); ++pid)
                {
                    auto const& old_cell = cell_particles[pid].iCell();
                    auto const pt        = tracker(old_cell, pps.local_tile_cell(old_cell));
                    per_particle(cell_particles[pid], layout, pt, pid, em);
                }
            }
        }
        else
        {
            auto const each = pps()[tile_idx]().size() / ws;

            auto const one = [&](std::size_t const pidx) {
                per_any_particle(parts, layout, tracker(parts.iCell(pidx), tile_cell), pidx, em);
            };

            for (std::size_t pid = 0; pid < each; ++pid)
                one(pid * ws + tidx);
        }
    }

    void per_any_particle(auto& particles, auto&&... args)
    {
        auto const& pidx = std::get<2>(std::forward_as_tuple(args...));
        per_particle(particles[pidx], args...);
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
        {
            using enum LayoutMode;
            if constexpr (Particles_t::layout_mode == AoSPCTS)
            {
                // AoSPCTS tracks per-cell buckets even in the tile ghost layer, so any
                // particle whose cell changed needs registering — domain and level ghost
                // alike. pt was built in per_tile, before advance() ran.
                pps.template move_check<particle_type>(pt, pidx, particle);
            }
            else if constexpr (particle_type == ParticleType::Domain)
            {
                if (isIn(particle, pps.box()))
                    pps.template move_check<particle_type>(pt, pidx, particle);
            }
        }
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
            else
            {
                using enum LayoutMode;
                // left the ghost box on the first half-step — register the departure or
                // the bucket holding it goes stale; move_check resolves it as a deletion
                if constexpr (boris_mode == MultiBorisMode::REF
                              and Particles_t::layout_mode == AoSPCTS)
                    pps.template move_check<particle_type>(pt, pidx, particle);
            }
        }
        else
        {
            PHARE_ASSERT(false);
        }
    }


    // one tile's worth of the COPY-mode deposit pass (AoSPCTS). Domain and LevelGhost
    // write into the SAME tile's rhoP/rhoC/flux, so callers must run this tile's Domain
    // and LevelGhost passes in one task (never two concurrent tasks for the same tile_idx)
    // — see MultiBorisPusherImplBase::move_cpu_copy_pooled. Different tiles are independent.
    void one_copy_tile_pcell(std::size_t const tile_idx, auto& boxings, auto& view, auto& pop)
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
                    if (isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                        interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                }
            }
        }
    }

    void per_copy_of_cpu_tile_pcell(auto& boxings, auto& view, auto& pop)
    {
        for (std::size_t tile_idx = 0; tile_idx < pps().size(); ++tile_idx)
            one_copy_tile_pcell(tile_idx, boxings, view, pop);
    }

    auto em_tile(auto const tidx)
    {
        return electromag.template as<Electromag_vt>([&](auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx](); });
        });
    }

    // one tile's worth of the COPY-mode deposit pass (AoSTS) — see one_copy_tile_pcell for
    // why Domain+LevelGhost of the SAME tile must stay in one task.
    void one_copy_tile(std::size_t const tile_idx, auto& boxings, auto& view, auto& pop)
    {
        auto const& patch_id = view.patchID();
        assert(boxings.count(patch_id));
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = pps.box();
        auto const is_border_     = [&](auto const& tile) {
            auto const tile_grow_box = grow(tile, 1);
            return patch_box * tile_grow_box != tile_grow_box;
        };

        auto& tile           = pps()[tile_idx];
        auto const tile_cell = pps.local_cell(tile.lower);
        auto const& layout   = electromag.E[0][tile_idx].layout();
        auto const em        = em_tile(tile_idx);
        auto& rhoP           = pop.particleDensity()[tile_idx];
        auto& rhoC           = pop.chargeDensity()[tile_idx];
        auto F = pop.flux().template as<VecField_vt>([&](auto& c) { return c()[tile_idx](); });
        Interpolator_t interp;

        auto on_tile = [&]<bool border>() {
            for (auto p : tile())
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

    void per_copy_of_cpu_tile(auto& boxings, auto& view, auto& pop)
    {
        for (std::size_t tile_idx = 0; tile_idx < pps().size(); ++tile_idx)
            one_copy_tile(tile_idx, boxings, view, pop);
    }

    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const dto2m;
    std::array<double, dim> halfdt;
};


// ── Common base: deduplicates move_rest / move_cpu_copy / sync_ref ───────────────────
// Derived must provide:
//   template<auto pt, auto mode> using Functors = MultiBorisFunctors<pt, mode, Derived>;
//   template<auto type> static void sync_particles(auto& particles);
//   using Particles_t = ...;

template<typename Derived>
struct MultiBorisPusherImplBase
{
    template<auto mode, typename ModelAccessor>
    static void move_rest(MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            typename Derived::template Functors<ParticleType::Domain, mode>{in, view, pop, domain,
                                                                            em}(in, i);

            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();
            typename Derived::template Functors<ParticleType::LevelGhost, mode>{
                in, view, pop, level_ghost, em}(in, i);
        }
    }

    template<auto mode, typename ModelAccessor>
    static void move_cpu_copy(MultiBoris<ModelAccessor>& in, auto& boxings, auto const i)
    {
        using enum LayoutMode;
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
            parts.reset_views();
            typename Derived::template Functors<particle_type, mode> fns{in, view, pop, parts, em};
            if constexpr (any_in(Derived::Particles_t::layout_mode, AoSPCTS))
                fns.per_copy_of_cpu_tile_pcell(boxings, view, pop);
            else
                fns.per_copy_of_cpu_tile(boxings, view, pop);
        };

        for (auto& pop : ions)
        {
            per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
            per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
        }
    }

    // ── Thread-pooled variants: level 1 (patch -> pool) is dispatched by the caller
    // (MultiBorisPusherImpl::move_pooled), which hands us the specific pool already
    // assigned to patch `i`. Here we do level 2: one task per tile, submitted to that
    // SAME pool so all of its threads_per_pool threads work this patch's tiles at once.
    // These must be called from the dispatching thread itself, never from inside a task
    // already running on `pool` — detach_task-then-wait on one's own pool deadlocks.

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

            using DomainFn = typename Derived::template Functors<ParticleType::Domain, mode>;
            using GhostFn  = typename Derived::template Functors<ParticleType::LevelGhost, mode>;

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
        using enum LayoutMode;
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();

            using DomainFn = typename Derived::template Functors<ParticleType::Domain, mode>;
            using GhostFn  = typename Derived::template Functors<ParticleType::LevelGhost, mode>;

            auto domain_fn = std::make_shared<DomainFn>(in, view, pop, domain, em);
            auto ghost_fn  = std::make_shared<GhostFn>(in, view, pop, level_ghost, em);

            auto const n_tiles = domain_fn->pps().size();
            assert(n_tiles == ghost_fn->pps().size());

            // view/pop are cheap resource-manager handles (ModelLevelAccessor::operator[]
            // already freezes an independent, patch-bound copy before returning -- see
            // amr/physical_models/models.hpp -- and pop's Field-type members are thin
            // repointable views), so capturing copies here is safe and correctly aliases
            // the real per-patch storage. Domain and LevelGhost for the SAME tile run
            // in the SAME task -- they write the same tile's rhoP/rhoC/flux, so they must
            // never be two concurrently-running tasks; different tiles are independent.
            for (std::size_t t = 0; t < n_tiles; ++t)
                pool.detach_task([domain_fn, ghost_fn, t, view, pop, &boxings]() mutable {
                    if constexpr (any_in(Derived::Particles_t::layout_mode, AoSPCTS))
                    {
                        domain_fn->one_copy_tile_pcell(t, boxings, view, pop);
                        ghost_fn->one_copy_tile_pcell(t, boxings, view, pop);
                    }
                    else
                    {
                        domain_fn->one_copy_tile(t, boxings, view, pop);
                        ghost_fn->one_copy_tile(t, boxings, view, pop);
                    }
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
            Derived::template sync_particles<ParticleType::Domain>(domain);

            auto& level_ghost = pop.levelGhostParticles();
            Derived::template sync_particles<ParticleType::LevelGhost>(level_ghost);
        }
    }
};


// ── MultiBorisPusherImpl specialisations for AoSTS ────────────────────────────────────

template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImplBase<MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles,
                                                    Electromag, Interpolator>>
{
    using This  = MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                                       Interpolator>;
    using Super = MultiBorisPusherImplBase<This>;

public:
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;


    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;

    template<auto type>
    static void sync_particles(auto& particles)
    {
        particles.template on_moved<type>();
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
                Super::template move_cpu_copy<mode>(in, boxings, i);
            else
                Super::template move_rest<mode>(in, i);

            if constexpr (not copy)
                Super::sync_ref(in, i);
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
                Super::template move_cpu_copy_pooled<mode>(pool, in, boxings, i);
            else
                Super::template move_rest_pooled<mode>(pool, in, i);
        }
        TP.sync(); // every patch's tiles, across every pool, are done

        if constexpr (not copy)
        {
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
                TP.async([&in, i] { Super::sync_ref(in, i); });
            TP.sync();
        }
    }
};


// AoSMapped: reuse AoSTS tile-based functors
template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSMapped, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                           Interpolator>
{
};


// AoSPCTS: tiles of per-cell AoS particles, tiled fields (same as AoSTS)
// Particle iteration differs: each tile contains a PerCellParticles container
// with flat iterators over all particles in that tile.
template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImplBase<MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout,
                                                    Particles, Electromag, Interpolator>>
{
    using This = MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout, Particles, Electromag,
                                      Interpolator>;
    using Super = MultiBorisPusherImplBase<This>;

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
                Super::template move_cpu_copy<mode>(in, boxings, i);
            else
                Super::template move_rest<mode>(in, i);

            if constexpr (not copy)
                Super::sync_ref(in, i);
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
                Super::template move_cpu_copy_pooled<mode>(pool, in, boxings, i);
            else
                Super::template move_rest_pooled<mode>(pool, in, i);
        }
        TP.sync(); // every patch's tiles, across every pool, are done

        if constexpr (not copy)
        {
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
                TP.async([&in, i] { Super::sync_ref(in, i); });
            TP.sync();
        }
    }
};

} // namespace PHARE::core


#endif
