// IWYU pragma: private, include "core/numerics/pusher/multi_boris.hpp"

#ifndef PHARE_CORE_PUSHER_BORIS_MKN_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_MULTI_BORIS_HPP

#include "core/data/particles/particle_array_def.hpp"


namespace PHARE::core
{

enum class MultiBorisMode : std::uint16_t { REF = 0, COPY };

// Primary (incomplete) template — specialisations provide the implementation per
// (LayoutMode, AllocatorMode) pair.  Add a new #include below for each new layout.
template<LayoutMode layout, AllocatorMode alloc, typename GridLayout, typename Particles,
         typename Electromag, typename Interpolator>
struct MultiBorisPusherImpl;


// ── MultiBoris state struct for AoSTS / ModelAccessor ─────────────────────────────────

template<typename ModelAccessor>
struct MultiBoris
{
    using Model_t      = ModelAccessor::Model_t;
    using GridLayout_t = ModelAccessor::GridLayout_t;

    static constexpr auto dim = GridLayout_t::dimension;
    using ParticleArray_t     = Model_t::particle_array_type;
    using Electromag_t        = Model_t::electromag_type;
    using Vecfield_t          = Electromag_t::vecfield_type;
    using Field_t             = Vecfield_t::field_type;
    using ParticleArray_v     = ParticleArray_t::view_t;
    using Box_t               = Box<int, dim>;
    using Particles_ptrs      = std::vector<ParticleArray_t*>;
    using StreamLauncher      = gpu::ThreadedStreamLauncher<ModelAccessor>;
    using PerTileParticles_t  = ParticleArray_t::per_tile_particles;

    MultiBoris(double const dt_, ModelAccessor& _accessor, std::function<void(int)> fn_ = {})
        : dt{dt_}
        , accessor{_accessor}
        , fn{fn_}
    {
    }

    void reset() {}

#if PHARE_HAVE_MKN_GPU
    using GpuBoxAlloc_t = mkn::gpu::ManagedAllocator<Box_t>;
#else
    using GpuBoxAlloc_t = std::allocator<Box_t>;
#endif
    using GpuBoxSpanSet_t = SpanSet<Box_t, default_span_size_t, GpuBoxAlloc_t>;

    double const dt;
    ModelAccessor& accessor;
    std::function<void(int)> fn;
    StreamLauncher streamer{accessor, 1};
    GpuBoxSpanSet_t gpu_nlgb;

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

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    static constexpr auto dim = GridLayout_t::dimension;
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    MultiBorisFunctors(auto& in, auto& view, auto& pop, auto& parts, auto& em)
        : pps{*parts}
        , electromag{em}
        , dto2m{0.5 * in.dt / pop.mass()}
        , halfdt{in.mesh(view.layout.meshSize(), in.dt)}
    {
        check_particles(parts);
        check_particles_views(parts);
    }

    void check(auto const& particle) _PHARE_ALL_FN_ { check_particle(particle); }

    void operator()(auto& in, [[maybe_unused]] auto const i)
    {
        if constexpr (Particles_t::alloc_mode == AllocatorMode::CPU)
            on_cpu_tiles(in);
        else if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_gpu_tiles(in, i);
        else
            throw std::runtime_error("Unknown alloc mode");
    }

    void on_gpu_tiles(auto& in, auto const i)
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
#endif
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
        auto const& layout                   = electromag.E[0][tile_idx].layout();
        auto const& em                       = em_tile(tile_idx);

        using enum LayoutMode;
        if constexpr (boris_mode == MultiBorisMode::REF and any_in(Particles_t::layout_mode))
        {
            aos_cpu_mapped_ref(tile, layout, em);
        }
        else
        {
            auto& parts          = tile();
            auto const each      = pps()[tile_idx]().size() / ws;
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

        check_electromag(em);
        check(particle);
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

    auto em_tile(auto const tidx) _PHARE_ALL_FN_
    {
        return electromag.template as<Electromag_vt>([&] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx](); });
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
    }

    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const dto2m;
    std::array<double, dim> halfdt;
};


// ── MultiBorisPusherImpl specialisations for AoSTS ────────────────────────────────────

template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
{
public:
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;
    using This = MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                                      Interpolator>;

    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;


    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelAccessor>
    static void move(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        static constexpr auto copy   = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu = Particles_t::alloc_mode == AllocatorMode::CPU;
        static constexpr auto is_gpu = Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED;

        if constexpr (copy and is_gpu)
        {
            std::vector<default_span_size_t> sizes;
            sizes.reserve(in.accessor.size());
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
                sizes.push_back(boxings.at(in.accessor[i].patchID()).nonLevelGhostBox.size());
            in.gpu_nlgb = typename MultiBoris<ModelAccessor>::GpuBoxSpanSet_t{std::move(sizes)};
            std::size_t offset = 0;
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
            {
                auto const& nlgb = boxings.at(in.accessor[i].patchID()).nonLevelGhostBox;
                std::copy(nlgb.begin(), nlgb.end(), in.gpu_nlgb.vec.begin() + offset);
                offset += nlgb.size();
            }
        }

        in.streamer.host([&](auto const i) mutable {
            if constexpr (copy and is_cpu)
                move_cpu_copy<mode>(in, boxings, i);
            else if constexpr (copy and is_gpu)
                move_gpu_copy<mode>(in, i);
            else
                move_rest<mode>(in, i);
        });

        in.streamer.host([&](auto const i) mutable {
            if constexpr (not copy)
                sync_ref(in, i);
        });
    }


    template<auto mode, typename ModelAccessor>
    static void move_cpu_copy(MultiBoris<ModelAccessor>& in, auto& boxings, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
            parts.reset_views();
            Functors<particle_type, mode> fns{in, view, pop, parts, em};
            fns.per_copy_of_cpu_tile(boxings, view, pop, parts);
        };

        for (auto& pop : ions)
        {
            per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
            per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
        }
    }


    template<auto mode, typename ModelAccessor>
    static void move_gpu_copy(MultiBoris<ModelAccessor>& in, auto const i)
    {
#if PHARE_HAVE_GPU
        using Tile_vt_ = Electromag_t::vecfield_type::field_type::value_type;
        using Interpolating_
            = Interpolating<Particles_t, GridLayout_t::interp_order,
                            (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)>;

        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto const dto2m      = 0.5 * in.dt / pop.mass();
            auto const halfdt     = MultiBoris<ModelAccessor>::mesh(view.layout.meshSize(), in.dt);
            auto const filter_box = in.gpu_nlgb[i];
            auto rhop             = pop.particleDensity();
            auto rhoc             = pop.chargeDensity();
            auto flux             = *pop.flux();
            auto const ds = static_cast<std::uint32_t>(ions.chargeDensity().max_tile_size());

            auto const launch = [&](auto parts) {
                if (parts().size() == 0)
                    return;
                using Launcher = gpu::ChunkLauncher<false>;
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = parts().size();
                launcher.ds
                    = ds * 5 * 8
                      + 2 * static_cast<std::uint32_t>(kernel::warp_size())
                            * static_cast<std::uint32_t>(sizeof(typename Particles_t::Particle_t))
                      + static_cast<std::uint32_t>(sizeof(int));
                assert(launcher.ds < 65000);
                launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                    Interpolating_::template on_tiles_push_deposit<Tile_vt_, ParticleType::Domain,
                                                                   Particles_t::alloc_mode>(
                        parts, em, flux, rhop, rhoc, filter_box, dto2m, halfdt);
                });
            };

            launch(*pop.domainParticles());
            launch(*pop.levelGhostParticles());
        }
#endif
    }


    template<auto mode, typename ModelAccessor>
    static void move_rest(MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (std::size_t j = 0; j < ions.size(); ++j)
        {
            auto& pop = ions[j];

            auto& domain = pop.domainParticles();
            domain.reset_views();
            Functors<ParticleType::Domain, mode>{in, view, pop, domain, em}(in, i);

            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();
            Functors<ParticleType::LevelGhost, mode>{in, view, pop, level_ghost, em}(in, i);
        }
    }


    template<typename ModelAccessor>
    static void sync_ref(MultiBoris<ModelAccessor>& in, auto const i)
    {
        if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            in.streamer.streams[i].sync();

        auto view      = in.accessor[i];
        auto [ions, _] = view.args;

        for (std::size_t j = 0; j < ions.size(); ++j)
        {
            auto& pop = ions[j];

            auto& domain = pop.domainParticles();
            sync_ts<ParticleType::Domain>(domain, in.streamer.streams[i]);

            auto& level_ghost = pop.levelGhostParticles();
            sync_ts<ParticleType::LevelGhost>(level_ghost, in.streamer.streams[i]);
        }
    }
};


// AoSMapped and AoSPC use the same tile-based functors as AoSTS
template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSMapped, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                           Interpolator>
{
};

template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSPC, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                           Interpolator>
{
};

} // namespace PHARE::core


#endif
