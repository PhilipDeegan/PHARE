#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/utilities/kernels.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"
#include "core/numerics/pusher/boris/basics.hpp"

#include "interpolator.hpp"

#include <string>
#include <iostream>
#include <optional>
#include <stdexcept>


namespace PHARE::core
{

// visits, for tile tidx, the particles of every cell its level ghost halo reaches, read
// from that cell's clamp-owner tile (TileSet::tag_cells_). Level ghosts are never in the
// patch domain, so only halo cells outside of it are visited: none for interior tiles
template<typename Particles>
void on_reachable_level_ghosts(Particles const& particles, std::size_t const tidx, auto&& fn)
    requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
{
    auto const& patch_box    = particles.box();
    auto const ghost_cells   = patch_box.lower[0] - particles.ghost_box().lower[0];
    auto const& tile         = particles()[tidx];
    auto const tile_halo     = grow(*tile, ghost_cells);
    auto const halo_in_patch = tile_halo * patch_box;
    if (halo_in_patch and *halo_in_patch == tile_halo)
        return; // interior tile

    for_each_point(tile_halo, [&](auto const& amr) {
        if (isIn(amr, patch_box))
            return;
        auto const& owner = (*particles().at(amr))();
        if constexpr (Particles::layout_mode == LayoutMode::AoSCMTS)
            for (auto const idx : owner.map(amr))
                fn(makeRange(owner, idx, idx + 1));
        else
            fn(owner(owner.local_cell(amr)));
    });
}


// how a per-tile GPU deposit kernel accumulates rhoP + the 3 flux components, all with atomics:
//  shared_all  all four tile fields staged in dynamic shared memory, one pass over the particles
//  shared_one  one field staged at a time, four passes (a pushing kernel pushes again per pass)
//  global      no staging, straight into the tile's fields in global memory
enum class TileDeposit : std::uint8_t { shared_all = 0, shared_one, global };

std::size_t constexpr tile_deposit_fields(TileDeposit const mode) _PHARE_ALL_FN_
{
    return mode == TileDeposit::shared_all ? 4 : mode == TileDeposit::shared_one ? 1 : 0;
}

// dynamic shared memory of a deposit kernel: staged fields of n_points each, then extra_bytes
std::size_t constexpr tile_deposit_bytes(TileDeposit const mode, std::size_t const n_points,
                                         std::size_t const extra_bytes = 0) _PHARE_ALL_FN_
{
    return tile_deposit_fields(mode) * n_points * sizeof(double) + extra_bytes;
}

// the fewest passes that fit in max_bytes of dynamic shared memory. For benchmarking,
//  PHARE_GPU_TILE_DEPOSIT=shared_all|shared_one|global forces a mode (throws if it can't fit)
inline TileDeposit tile_deposit_mode(std::size_t const n_points, std::size_t const extra_bytes,
                                     std::size_t const max_bytes)
{
    using enum TileDeposit;
    static std::optional<TileDeposit> const forced = []() -> std::optional<TileDeposit> {
        auto const env = core::get_env("PHARE_GPU_TILE_DEPOSIT");
        if (!env)
            return std::nullopt;
        if (*env == "shared_all")
            return shared_all;
        if (*env == "shared_one")
            return shared_one;
        if (*env == "global")
            return global;
        throw std::runtime_error("PHARE_GPU_TILE_DEPOSIT: unknown mode " + *env);
    }();

    auto const fits = [&](TileDeposit const mode) {
        return tile_deposit_bytes(mode, n_points, extra_bytes) <= max_bytes;
    };

    if (forced)
    {
        if (!fits(*forced))
        {
            // printed too: this runs on stream worker threads, where the exception terminates
            std::string const msg
                = "PHARE_GPU_TILE_DEPOSIT: forced mode needs "
                  + std::to_string(tile_deposit_bytes(*forced, n_points, extra_bytes))
                  + " bytes of shared memory, the maximum is " + std::to_string(max_bytes);
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }
        return *forced;
    }
    for (auto const mode : {shared_all, shared_one})
        if (fits(mode))
            return mode;
    return global;
}

// calls fn(std::integral_constant<TileDeposit, mode>{}), to make a runtime mode a kernel template
//  argument
template<typename Fn>
void with_tile_deposit(TileDeposit const mode, Fn&& fn)
{
    using enum TileDeposit;
    if (mode == shared_all)
        fn(std::integral_constant<TileDeposit, shared_all>{});
    else if (mode == shared_one)
        fn(std::integral_constant<TileDeposit, shared_one>{});
    else
        fn(std::integral_constant<TileDeposit, global>{});
}


// simple facade to launch e.g. GPU kernels if needed
//  and to not need to modify the interpolator much for gpu specifically

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false,
         typename Interpolator_t = Interpolator<dim, interpOrder, atomic_ops>>
class Interpolating
{
public:
    template<typename Particles>
    inline void operator()(Particles const& particles, auto& rhoP, auto& flux, auto const& layout,
                           double coef = 1.)
    {
        PHARE_LOG_SCOPE(3, "Interpolating::particleToMesh");
        particleToMesh(particles, layout, rhoP, flux, coef);
    }

    // LevelGhost: level ghosts live only in their clamp-owner tile and are deposited before
    // the tile reduction + border sums (fillFluxBorders/fillDensityBorders), so like domain
    // particles each tile deposits only its own, the reduction completes the tile seams
    // GPU_UNIFIED falls back to the CPU path until there's an optimized GPU version
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    _PHARE_ALL_FN_
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
            particleToMesh_tile<type>(tidx, particles, rhoP, flux, coef);
    }

    // only writes the moments of tile tidx, distinct tiles can be deposited concurrently
    // AoSTS GPU_UNIFIED: host loop over unified memory like its particleToMesh below, the
    //  device deposit is the updater's (on_tiles / on_tiles_push_deposit)
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh_tile(std::size_t const tidx, Particles const& particles, auto& rhoP,
                             auto& flux, double coef = 1.)
        requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS,
                        LayoutMode::AoSTS))
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        auto [rhop, F]      = tiles_at(tidx, rhoP, flux);
        auto const& rho_lay = tile_layout(rhoP, tidx);
        auto const deposit  = [&](auto const& cell_particles) {
            for (auto const& p : cell_particles)
                interp_.particleToMesh(p, rhop, F, rho_lay, coef);
        };

        auto& tile_particles = particles()[tidx]();
        if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSCMTS, LayoutMode::AoSTS))
            deposit(tile_particles);
        else
            for_each_point(tile_particles.local_box(),
                           [&](auto const& bix) { deposit(tile_particles(bix)); });
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSMapped)
    _PHARE_ALL_FN_
    {
        interp_(particles, rhoP, flux, layout, coef);
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSTS)
    _PHARE_ALL_FN_
    {
        // GPU_UNIFIED falls back to the CPU path until there's an optimized GPU version
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto [rhop, F]      = tiles_at(tidx, rhoP, flux);
            auto const& rho_lay = tile_layout(rhoP, tidx);
            for (auto const& p : particles()[tidx]())
                interp_.particleToMesh(p, rhop, F, rho_lay, coef);
        }
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::SoATS)
    _PHARE_ALL_FN_
    {
        for (auto const& tile : particles())
            for (std::size_t i = 0; i < tile().size(); ++i)
                interp_.particleToMesh(tile()[i], rhoP, flux, layout, coef);
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSPC)
    _PHARE_ALL_FN_
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;

        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            for_each_point(particles.local_box(), [&](auto const& bix) {
                for (auto const& p : particles(bix))
                    interp_.particleToMesh(p, rhoP, flux, layout, coef);
            });
        }
        // else if (alloc_mode == AllocatorMode::GPU_UNIFIED and impl < 2)
        // {
        //     static_assert(atomic_ops, "GPU must be atomic");
        //     PHARE_WITH_MKN_GPU(
        //         mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
        //             auto it = particles[mkn::gpu::idx()];
        //             Interpolator_t{}.particleToMesh(*it, density, flux, layout, coef);
        //         }); //
        //     )
        // }
        else if (alloc_mode == AllocatorMode::GPU_UNIFIED /*and impl == 2*/)
        {
            PHARE_WITH_MKN_GPU({ //
                using lobox_t  = Particles::lobox_t;
                using Launcher = gpu::BoxCellNLauncher<lobox_t>;
                auto lobox     = particles.local_box();
                Launcher{lobox, particles.max_size()}([=] _PHARE_ALL_FN_() mutable {
                    box_kernel(particles, layout, flux, rhoP, coef);
                });
            })
        }
        else
            throw std::runtime_error("fail");
    }

    // stream: optional, GPU_UNIFIED SoAVX deposits on it (queued after prior work on it)
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1., auto&... stream)
        requires(any_in(Particles::layout_mode, LayoutMode::AoS, LayoutMode::SoA,
                        LayoutMode::SoAVX, LayoutMode::SoAVXChunked))
    _PHARE_ALL_FN_
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;

        // not atomic (level init, messenger ghosts): on the host, GPU_UNIFIED memory included
        if constexpr (alloc_mode == AllocatorMode::CPU or !atomic_ops)
        {
            if constexpr (is_soavx(Particles::layout_mode))
            {
                for (auto const& it : particles)
                    interp_.particleToMesh(it, rhoP, flux, layout, coef);
            }
            else
                interp_(particles, rhoP, flux, layout, coef);
        }
        else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            static_assert(atomic_ops, "GPU must be atomic");
            if constexpr (is_soavx(Particles::layout_mode))
            {
                // a thread per particle, loaded from the columns, global atomics
                PHARE_WITH_MKN_GPU( //
                    SoAVXConstSpanParticles<Particles::dimension, alloc_mode> const ps{
                        particles}; //
                    kernel::launch(stream..., ps.size(),
                                   [=] _PHARE_ALL_FN_() mutable {
                                       Interpolator_t{}.particleToMesh(ps.copy(kernel::idx()),
                                                                       rhoP, flux, layout, coef);
                                   }); //
                )
            }
            else
            {
                PHARE_WITH_MKN_GPU( //
                    kernel::launch(particles.size(),
                                   [=] _PHARE_ALL_FN_() mutable {
                                       auto particle = particles.begin() + kernel::idx();
                                       Interpolator_t{}.particleToMesh(deref(particle), rhoP, flux,
                                                                       layout, coef);
                                   }); //
                )
            }
        }
        else
            throw std::runtime_error("fail");
    }

    // AoSPCTS GPU deposit, launched with gpu::TileCellLauncher: block = (tile, cell), the block's
    // threads stride the cell's particles. Cells share field nodes, so the interpolator must be
    // atomic. Field_t as on_tiles (tiled field value_type)
    template<typename Field_t, typename Launcher>
    static void on_tile_cells(auto& particles, auto& Flux, auto& Rhop, std::size_t const threads,
                              double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        static_assert(atomic_ops, "GPU must be atomic");
        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;

        auto const tidx   = Launcher::tile_idx();
        auto& cps         = particles()[tidx]();
        auto const& lobox = cps.local_box();
        if (Launcher::cell_idx() >= lobox.size())
            return;

        auto& parts        = cps(point_at(lobox, Launcher::cell_idx()));
        auto& rhop_tile    = Rhop[tidx];
        auto flux          = Flux.template as<TensorField_t>([&](auto& c) { return c[tidx](); });
        auto const& layout = rhop_tile.layout();
        auto& pdensity     = rhop_tile();
        for (std::size_t pid = Launcher::thread_idx(); pid < parts.size(); pid += threads)
            Interpolator_t{}.particleToMesh(parts[pid], pdensity, flux, layout, coef);
#endif // PHARE_HAVE_MKN_GPU_HW
    }

    template<bool in_box = false, typename Particles, typename GridLayout, typename VecField,
             typename Field>
    static void box_kernel(Particles& particles, GridLayout& layout, VecField& flux, Field& density,
                           double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        using lobox_t         = Particles::lobox_t;
        using Launcher        = gpu::BoxCellNLauncher<lobox_t>;
        auto const& lobox     = particles.local_box();
        auto const& blockidx  = Launcher::block_idx();
        auto const& threadIdx = Launcher::thread_idx();
        auto& parts           = particles(point_at(lobox, blockidx));

        auto const interp = [&]() {
            Interpolator_t{}.particleToMesh(parts[threadIdx], density, flux, layout, coef);
        };

        if (threadIdx < parts.size())
        {
            if constexpr (!in_box)
            {
                interp();
            }
            else
            {
                if (isIn(parts[threadIdx], particles.box()))
                    interp();
            }
        }
#endif // PHARE_HAVE_MKN_GPU_HW
    }

#if 0 // old unused
    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel_ts(Particles& particles, GridLayout& layout, VecField& flux,
                                Field& density, double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW

        Point<std::uint32_t, 3> const rcell{threadIdx.z, threadIdx.y, threadIdx.x};
        extern __shared__ double data[];

        auto& tile        = *particles().at(rcell + 2);
        auto& parts       = tile();
        auto const each   = parts.size() / 8;
        auto const tile_p = rcell % 2;
        auto const pidx   = tile_p[2] + tile_p[1] * 2 + tile_p[0] * 2 * 2;
        auto const ji     = [&](auto const i) { return i * 8 + pidx; };
        auto doX          = [&]<std::uint8_t IDX = 0>(auto& feeld) {
            auto v = make_array_view(&data[0], feeld.shape());
            if (kernel::idx() == 0)
                for (auto& e : v)
                    e = 0;
            Interpolator_t interp;
            std::size_t i = 0, j = ji(i);
            __syncthreads();
            for (; i < each; ++i, j = ji(i))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            if (pidx < parts.size() - (8 * each))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            __syncthreads();
            if (kernel::idx() == 0)
                for (i = 0; i < v.size(); ++i)
                    feeld.data()[i] = v.data()[i];
        };

        doX(density);
        doX.template operator()<1>(flux[0]);
        doX.template operator()<2>(flux[1]);
        doX.template operator()<3>(flux[2]);


#endif // PHARE_HAVE_MKN_GPU_HW
    }
#endif

    // accumulates into one tile's rhoP + flux (global memory) as deposit_v says, see TileDeposit.
    //  each_particle(dep) must call dep(particle) for every particle this thread deposits; every
    //  thread of the block calls it, and it may __syncthreads. With shared_one it runs four times
    template<auto deposit_v>
    static void tile_deposit(auto pdensity, auto flux, auto const& tbox, double const coef,
                             auto&& each_particle) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        static_assert(atomic_ops, "GPU must be atomic");
        using enum TileDeposit;

        extern __shared__ double data[];
        auto const ziz  = static_cast<int>(product(tbox.shape()));
        auto const tidx = static_cast<int>(threadIdx.x);
        auto const bdim = static_cast<int>(blockDim.x);

        // pdensity/flux are copies of the tile's views (same data), rho/fs are restaged
        auto rho = pdensity;
        auto fs  = flux;

        // staged sums into the tile's fields: several blocks per tile (gridDim.y) add atomically
        auto const flush = [multi = gridDim.y > 1](double& dst, double const v) {
            if (multi)
                atomicAdd(&dst, v);
            else
                dst += v;
        };

        if constexpr (deposit_v == global)
        {
            each_particle([&](auto const& p) {
                Interpolator_t interp;
                interp.particleToMesh(p, rho, fs, tbox, coef);
            });
        }
        else if constexpr (deposit_v == shared_all)
        {
            auto s_rho = make_array_view(&data[ziz * 0], pdensity.shape());
            auto s_fx  = make_array_view(&data[ziz * 1], flux[0].shape());
            auto s_fy  = make_array_view(&data[ziz * 2], flux[1].shape());
            auto s_fz  = make_array_view(&data[ziz * 3], flux[2].shape());

            for (int i = tidx; i < ziz * 4; i += bdim)
                data[i] = 0;
            __syncthreads();

            rho.reset(s_rho);
            fs[0].reset(s_fx);
            fs[1].reset(s_fy);
            fs[2].reset(s_fz);
            each_particle([&](auto const& p) {
                Interpolator_t interp;
                interp.particleToMesh(p, rho, fs, tbox, coef);
            });
            __syncthreads();

            for (int i = tidx; i < ziz; i += bdim)
            {
                flush(pdensity.data()[i], s_rho.data()[i]);
                flush(flux[0].data()[i], s_fx.data()[i]);
                flush(flux[1].data()[i], s_fy.data()[i]);
                flush(flux[2].data()[i], s_fz.data()[i]);
            }
        }
        else if constexpr (deposit_v == shared_one)
        {
            for_N<4>([&](auto ic) {
                auto constexpr c = ic();
                auto& dst        = [&]() -> auto& {
                    if constexpr (c == 0)
                        return pdensity;
                    else
                        return flux[c - 1];
                }();
                auto staged = dst;
                auto s_f    = make_array_view(&data[0], dst.shape());

                for (int i = tidx; i < ziz; i += bdim)
                    data[i] = 0;
                __syncthreads();

                staged.reset(s_f);
                each_particle([&](auto const& p) {
                    Interpolator_t interp;
                    interp.p2m_setup(p, tbox);
                    interp.template p2m_per_component<c>(p, staged, coef);
                });
                __syncthreads();

                for (int i = tidx; i < ziz; i += bdim)
                    flush(dst.data()[i], s_f.data()[i]);
                __syncthreads(); // shared memory is reused by the next field
            });
        }
#endif
    }


    template<typename Field_t, auto deposit_v = TileDeposit::shared_all>
    static void on_tiles(auto& particles, auto& Flux, auto& Rhop, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;

        auto& rhop_tile = Rhop[blockIdx.x];
        auto const ft = Flux.template as<TensorField_t>([&](auto& c) { return c[blockIdx.x](); });
        auto& parts   = particles()[blockIdx.x]();

        tile_deposit<deposit_v>(rhop_tile(), ft, rhop_tile.ghost_box(), coef, [&](auto&& dep) {
            for (std::size_t i = threadIdx.x; i < parts.size(); i += blockDim.x)
                dep(parts[i]);
        });
#endif
    }

    // dynamic shared memory on_tiles_push_deposit needs after its staged fields, for block_size
    //  threads per block
    template<typename Particle_t>
    static std::size_t push_deposit_extra_bytes(std::size_t const block_size)
    {
        return 2 * block_size * sizeof(Particle_t) + sizeof(int);
    }

    // Combined Boris push + particle-to-mesh deposit in a single GPU kernel.
    // Originals are untouched: each thread works on a register copy.
    // Qualifying particles (isIn filter_box and not in exclude_boxes after push) are atomically
    // placed into a shared block buffer. Once the buffer holds a full block the whole block
    // deposits together; any remainder is flushed at the end. With shared_one the push runs
    // once per staged field.
    // Grid: x = tile, y = blocks per tile, each block a strided share of the tile's particles.
    //
    // Shared memory layout (launcher.ds = tile_deposit_bytes(deposit_v, ziz,
    //                                                        push_deposit_extra_bytes<Particle_t>())):
    //   [tile_deposit_fields(deposit_v) * ziz * sizeof(double)]  staged fields
    //   [2 * blockDim.x * sizeof(Particle_t)]                    block buffer
    //   [sizeof(int)]                                            atomic fill counter
    //
    // 2*blockDim.x slots are the minimum needed: before a batch buf_count can be
    // blockDim.x-1, and if every thread in the batch qualifies buf_count reaches
    // 2*blockDim.x-1 — still within bounds.
    template<typename Field_t, auto particle_type_v, auto alloc_mode_v,
             auto deposit_v = TileDeposit::shared_all>
    static void on_tiles_push_deposit(auto& particles, auto const& electromag, auto& Flux,
                                      auto& Rhop, auto const& filter_box, auto const& exclude_boxes,
                                      double const coef1, auto const& halfdt,
                                      double const coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;
        using Scalar_vt     = typename Field_t::value_type;
        using VecField_vt   = basic::TensorField<Scalar_vt, 1>;
        using Electromag_vt = basic::Electromag<VecField_vt>;
        using Particle_t    = ParticleDefaults<dim>::Particle_t;

        extern __shared__ double data[];

        auto& rhop_tile = Rhop[blockIdx.x];
        auto const ft   = Flux.template as<TensorField_t>([&](auto& c) { return c[blockIdx.x](); });
        auto const tbox = rhop_tile.ghost_box();
        auto const ziz  = product(tbox.shape());
        auto const tidx = threadIdx.x;
        auto const ws   = blockDim.x;

        // block buffer and atomic counter after the staged fields
        auto* const part_buf = reinterpret_cast<Particle_t*>(
            reinterpret_cast<char*>(&data[0])
            + tile_deposit_fields(deposit_v) * ziz * sizeof(double));
        auto& buf_count = *reinterpret_cast<int*>(reinterpret_cast<char*>(part_buf)
                                                  + 2 * ws * sizeof(Particle_t));

        // Per-tile EM view and layout
        auto const em = electromag.template as<Electromag_vt>([=] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>(
                [&vf] _PHARE_ALL_FN_(auto i) { return vf[i][blockIdx.x](); });
        });
        auto const& layout = electromag.E[0][blockIdx.x].layout();

        auto& parts        = particles()[blockIdx.x]();
        auto const n       = parts.size();
        auto const stride  = static_cast<std::size_t>(ws) * gridDim.y; // a batch over all blocks
        auto const nbatch  = (n + stride - 1) / stride;
        auto const b_start = static_cast<std::size_t>(blockIdx.y) * ws;

        auto const push_and_deposit = [&](auto&& dep) {
            if (tidx == 0)
                buf_count = 0;
            __syncthreads();

            for (std::size_t batch = 0; batch < nbatch; ++batch)
            {
                std::size_t const pidx = batch * stride + b_start + tidx;
                if (pidx < n)
                {
                    Particle_t p = parts[pidx]; // register copy — original untouched

                    p.iCell() = boris::advance<alloc_mode_v>(p, halfdt);
                    {
                        Interpolator_t interp;
                        boris::accelerate(p, interp.m2p(p, em, layout), coef1);
                    }
                    p.iCell() = boris::advance<alloc_mode_v>(p, halfdt);

                    if (isIn(p.iCell(), filter_box) and !isIn(p.iCell(), exclude_boxes))
                    {
                        int const slot = atomicAdd(&buf_count, 1);
                        part_buf[slot] = p;
                    }
                }

                __syncthreads();

                // Once a full warp is available, deposit it and compact the buffer.
                if (buf_count >= static_cast<int>(ws))
                {
                    dep(part_buf[tidx]);

                    __syncthreads();

                    int const remaining = buf_count - static_cast<int>(ws);
                    if (static_cast<int>(tidx) < remaining)
                        part_buf[tidx] = part_buf[tidx + ws];
                    if (tidx == 0)
                        buf_count = remaining;

                    __syncthreads();
                }
            }

            // Flush any remaining particles (partial warp).
            if (static_cast<int>(tidx) < buf_count)
                dep(part_buf[tidx]);
        };

        tile_deposit<deposit_v>(rhop_tile(), ft, tbox, coef, push_and_deposit);
#endif
    }


    template<typename Field, typename Particles, typename GridLayout>
    static void on_cpu_tiles(Particles& particles, GridLayout& layout,
                             double coef = 1.) _PHARE_ALL_FN_
    {
        using TensorField_t = basic::TensorField<typename Field::Super, 1>;

        for (auto& tile : particles())
        {
            auto& parts                                = tile();
            auto const tbox                            = tile.field_box();
            auto const& [density, flux0, flux1, flux2] = tile.fields();
            auto const r0                              = density;

            TensorField_t flux{flux0, flux1, flux2};
            Interpolator_t interp;
            std::size_t pid = 0;

            interp.particleToMesh(parts[pid], density, flux, /*layout,*/ tbox, coef);
        }
    }




    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel(Particles& particles, GridLayout& layout, VecField& flux,
                             Field& density, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        auto constexpr static nghosts = GridLayout::options.field_ghost_width;
        extern __shared__ double data[];
        auto const& lobox = particles.local_box();
        auto const ziz    = 9 * 9 * 9;

        auto const t_x = threadIdx.x;
        auto const t_y = threadIdx.y;
        auto const t_z = threadIdx.z;

        core::Point<std::uint32_t, 3> const tcell{t_x, t_y, t_z};
        core::Point<std::uint32_t, 3> const locell = tcell + nghosts;

        {
            Interpolator_t interp;
            auto& parts = particles(locell);

            { // rho
                auto const r0 = density.data();
                auto v        = make_array_view(&data[0], density.shape());
                density.setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], v);
                }
                __syncthreads();
                density.setData(r0);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        density.data()[i] = v.data()[i];
            }

            {
                auto const f0 = flux[0].data();
                auto v        = make_array_view(&data[0], flux[0].shape());
                flux[0].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[0]);
                }
                __syncthreads();
                flux[0].setData(f0);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[0].data()[i] = v.data()[i];
            }

            {
                auto const f1 = flux[1].data();
                auto v        = make_array_view(&data[0], flux[1].shape());
                flux[1].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[1]);
                }
                __syncthreads();
                flux[1].setData(f1);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[1].data()[i] = v.data()[i];
            }

            {
                auto const f2 = flux[2].data();
                auto v        = make_array_view(&data[0], flux[2].shape());
                flux[2].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[2]);
                }
                __syncthreads();
                flux[2].setData(f2);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[2].data()[i] = v.data()[i];
            }
        }


#endif // PHARE_HAVE_MKN_GPU_HW
    }




    Interpolator_t interp_;
};

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
struct MomentumTensorInterpolating
{
    using Interpolator_t = MomentumTensorInterpolator<dim, interpOrder, atomic_ops>;

    Interpolator_t interp_;

public:
    template<typename Particles_t>
    inline void operator()(Particles_t& particles, auto& momentumTensor, auto const& layout,
                           double mass = 1.)
    {
        assert(momentumTensor.isUsable());
        particleToMesh(particles, momentumTensor, layout, mass);
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSMapped)
    {
        interp_(particles, momentumTensor, layout, mass);
    }

    // SoAVX columns have no particle reference to iterate: deposited from an AoS buffer
    //  loaded in batches, on the host (diagnostics only)
    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(is_soavx(Particles_t::layout_mode))
    {
        static constexpr std::size_t batch = 1024;
        using Buffer_t                     = ParticleArray<ParticleArrayOptions{
                                Particles_t::dimension, LayoutMode::AoS}>;

        Buffer_t buffer;
        buffer.reserve(std::min(batch, particles.size()));
        for (std::size_t i = 0; i < particles.size(); ++i)
        {
            buffer.emplace_back(particles.copy(i));
            if (buffer.size() == batch or i + 1 == particles.size())
            {
                interp_(buffer, momentumTensor, layout, mass);
                buffer.clear();
            }
        }
    }

    // unlike Interpolating::particleToMesh, level ghosts are deposited after the tile reduction
    // (fillPopMomTensor), so each tile gathers every level ghost its halo reaches
    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);
            auto const deposit
                = [&](auto const& cell_particles) { interp_(cell_particles, mt, mt_lay, mass); };

            if constexpr (type == ParticleType::LevelGhost)
                on_reachable_level_ghosts(particles, tidx, deposit);
            else
            {
                auto& cps = particles()[tidx]();
                for_each_point(cps.local_box(), [&](auto const& bix) { deposit(cps(bix)); });
            }
        }
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSCMTS)
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);
            auto const deposit
                = [&](auto const& tile_particles) { interp_(tile_particles, mt, mt_lay, mass); };

            if constexpr (type == ParticleType::LevelGhost)
                on_reachable_level_ghosts(particles, tidx, deposit);
            else
                deposit(particles()[tidx]());
        }
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSTS)
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt = tile_at(momentumTensor, tidx); // interp_ takes it by reference
            interp_(particles()[tidx](), mt, tile_layout(momentumTensor, tidx), mass);
        }
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
