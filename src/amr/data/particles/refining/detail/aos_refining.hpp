// IWYU pragma: private, include "amr/data/particles/refining/particles_refining.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER


#include "core/data/particles/particle_array_def.hpp"
#include "core/def.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_selector.hpp"
#include "core/operators.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box_spans.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/refining/detail/def_refining.hpp"

#include <SAMRAI/hier/BoxContainer.h>

#include <vector>


namespace PHARE::amr
{

using enum core::LayoutMode;
using enum AllocatorMode;


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    reserve_flat<type, options>(src, dst, dst_boxes, dst_amr_box);
    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box = phare_box_from<Src::dimension>(samrai_box);
        stream_split_particles<type, options>(src, dst, dst_box, fn0, fn1);
    }
}

template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSMapped, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;

    auto const& dstBox   = dst.box();
    auto const ghost_nbr = (dst.ghost_box().shape()[0] - dst.box().shape()[0]) / 2;

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const& splitBox    = dst_box;
        auto const coarseDstBox = coarsen_box(dst_box);

        for (auto const& src_tile : src())
        {
            if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                continue;

            if constexpr (type == core::ParticleType::Domain)
            {
                for (auto& dst_tile : dst())
                    if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                        if (auto const overlap = dst_box * dst_tile)
                            reserve_flat<type, options>(src_tile(), dst_tile(), *overlap,
                                                        dst_amr_box, 1.);

                // DOUBLE LOOP == LESS COSTLY LONG TERM
                for (auto& dst_tile : dst())
                    if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                        if (auto const overlap = dst_box * dst_tile)
                            stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap,
                                                                  fn0, fn1);
            }
            else if constexpr (type == core::ParticleType::Ghost)
            {
                for (auto& dst_tile : dst())
                {
                    // Only grow in patch-boundary directions to avoid double-counting ghost
                    // cells near tile boundaries in 2D/3D (adjacent tiles share grow overlap).
                    auto excl_box = grow(dst_tile, ghost_nbr);
                    for (std::size_t d = 0; d < Src::dimension; ++d)
                    {
                        if (dst_tile.lower[d] != dstBox.lower[d])
                            excl_box.lower[d] = dst_tile.lower[d];
                        if (dst_tile.upper[d] != dstBox.upper[d])
                            excl_box.upper[d] = dst_tile.upper[d];
                    }
                    if (*(dst.box() * excl_box) != excl_box)
                        if (auto const overlap = dst_box * excl_box)
                        {
                            reserve_flat<type, options>(src_tile(), dst_tile(), *overlap,
                                                        dst_amr_box, 1.);
                            stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap,
                                                                  fn0, fn1);
                        }
                }
            }
            else
            {
                assert(false);
            }
        }
    }
}



template<auto type, typename Dst, typename Tile, typename Box_t>
std::optional<Box_t> refine_tile_overlap(Dst const& dst, Tile const& dst_tile, Box_t const& dst_box,
                                         std::size_t const ghost_nbr)
{
    // dst_box can extend past this array's own domain (periodic wrap, a restricted/
    // temporary ParticlesData SAMRAI builds internally, or the level ghost layer).
    // Particles are never duplicated across tiles, so each bound's margin is owned
    // by whichever tile touches the domain edge it's adjacent to (same clamp rule
    // tag_cells_ uses, so TileSet::at(ghost_cell) finds them); a tile that doesn't
    // own a bound just clips to its own edge there instead of the whole region
    // being rejected. Readers needing a neighbour's level ghosts gather via at().
    auto const growbox = grow(dst_tile, ghost_nbr);
    if (not(dst_box * growbox))
        return std::nullopt;

    auto const& domain_box = dst.box();
    auto const& tile_box   = *dst_tile;
    auto region            = dst_box;
    for (std::size_t d = 0; d < Dst::dimension; ++d)
    {
        if (dst_box.lower[d] < domain_box.lower[d])
        {
            bool const owns_edge = domain_box.lower[d] >= tile_box.lower[d]
                                   and domain_box.lower[d] <= tile_box.upper[d];
            region.lower[d]
                = owns_edge ? dst_box.lower[d] : std::max(dst_box.lower[d], tile_box.lower[d]);
        }
        else
            region.lower[d] = std::max(dst_box.lower[d], tile_box.lower[d]);

        if (dst_box.upper[d] > domain_box.upper[d])
        {
            bool const owns_edge = domain_box.upper[d] >= tile_box.lower[d]
                                   and domain_box.upper[d] <= tile_box.upper[d];
            region.upper[d]
                = owns_edge ? dst_box.upper[d] : std::min(dst_box.upper[d], tile_box.upper[d]);
        }
        else
            region.upper[d] = std::min(dst_box.upper[d], tile_box.upper[d]);

        if (region.lower[d] > region.upper[d])
            return std::nullopt;
    }
    return region;
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSPCTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // per-tile destination is per-cell storage, with no single flat buffer for reserve()
    // to size - see the shrink_to_fit_in comment below.
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    auto const ghost_nbr                      = options.particle_ghost_width;

    // Per-cell reserve: each dst cell's own coarse-parent cell count (not a box-wide average)
    // - uneven patterns (3D's 6-particle Pink) still spread unevenly across a coarse cell's
    // children, but that's a per-cell rounding slop, not the systematic box-total error a
    // uniform box-average would introduce.
    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const coarseDstBox = coarsen_box(dst_box);

        auto const dst_tile_overlap = [&](auto& dst_tile) {
            return refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr);
        };

        // reserve pass: one .reserve() call per dst cell, via reserve_flat_count against the
        // whole (tiled) src - count_particles<AoSPCTS> sums each tile's own domain box only
        // (never ghost cells), so cells sitting in more than one neighbor tile's ghost halo
        // don't get summed twice the way checking every src_tile's ghost_box used to.
        for (auto& dst_tile : dst())
            if (auto const overlap = dst_tile_overlap(dst_tile))
            {
                auto& dst_pc = dst_tile();
                core::for_each_point(*overlap, [&](auto const& bix) {
                    core::Box<int, Src::dimension> const cell_box{*bix, *bix};
                    dst_pc(dst_pc.local_cell(*bix))
                        .reserve(reserve_flat_count<type, options>(
                            src, cell_box, dst_amr_box,
                            per_cell_density_factor<options, type, Src>()));
                });
            }

        for (auto const& src_tile : src())
        {
            if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                continue;

            for (auto& dst_tile : dst())
                if (auto const overlap = dst_tile_overlap(dst_tile))
                    stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap, fn0,
                                                          fn1);
        }
    }

    dst.template on_appended<type>();
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSCMTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // one flat cell-mapped array per tile: each dst tile has a single flat buffer, so
    // reserve per tile overlap
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    auto const ghost_nbr                      = options.particle_ghost_width;

    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    using Box_t     = core::Box<int, Src::dimension>;
    using DstTile_t = std::decay_t<decltype(dst()[0])>;
    struct DstOverlap
    {
        DstTile_t* tile;
        Box_t overlap;      // fine cells of dst_box this dst tile receives
        Box_t coarse_reach; // coarse cells whose particles can split into overlap
    };
    thread_local std::vector<DstOverlap> dst_overlaps; // reused, no per call allocation

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const coarseDstBox = coarsen_box(dst_box);

        dst_overlaps.clear();
        for (auto& dst_tile : dst())
            if (auto const overlap = refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr))
            {
                reserve_flat<type, options>(src, dst_tile(), *overlap, dst_amr_box,
                                            per_cell_density_factor<options, type, Src>());
                // same coarse box stream_split_particles filters source particles with
                dst_overlaps.emplace_back(&dst_tile, *overlap,
                                          coarsen_box(grow(*overlap, ghost_nbr)));
            }

        for (auto const& src_tile : src())
        {
            auto const src_reach = grow(src_tile, ghost_nbr); // cells src_tile() may hold
            if (not(coarseDstBox * src_reach))
                continue;

            // only dst tiles whose coarse reach touches src_tile can receive its children,
            // skip scanning src_tile's particles for every other dst tile
            for (auto const& [dst_tile, overlap, coarse_reach] : dst_overlaps)
                if (coarse_reach * src_reach)
                    stream_split_particles<type, options>(src_tile(), (*dst_tile)(), overlap, fn0,
                                                          fn1);
        }
    }

    dst.template on_appended<type>();
}



template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSMapped, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                          auto fn1)
{
    ParticlesRefiner<AoSMapped, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


// one block per (dst box, src tile reaching it): each thread splits a coarse particle, fn1 keeps
//  its children in the dst box (see device_refining), gathered in managed memory sized for every
//  child of every listed particle, then appended to the dst tiles on the device
template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSTS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                      auto fn1)
{
    static_assert(device_refining<Src>);
    PHARE_LOG_SCOPE(3, "ParticlesRefiner<AoSTS, GPU_UNIFIED>::operator()");

#if PHARE_HAVE_MKN_GPU
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;

    auto constexpr dim    = Src::dimension;
    auto constexpr n_kids = decltype(fn1)::max_children;
    using Box_t           = core::Box<int, dim>;
    using Particle_t      = typename Src::Particle_t;
    using Launcher        = core::gpu::BlockPerLauncher<true>;
    using Op              = core::Operators<unsigned long long, true, true>;

    struct TileSplit
    {
        std::size_t tile;
        Box_t dst_box;
        Box_t split_box;  // fine cells a child of dst_box can come from
        Box_t coarse_box; // its coarse cells
    };

    auto const ghost_nbr = options.particle_ghost_width;
    std::vector<TileSplit> host_splits;
    std::size_t n_src = 0;
    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box    = phare_box_from<dim>(samrai_box);
        auto const split_box  = grow(dst_box, ghost_nbr);
        auto const coarse_box = coarsen_box(split_box);
        for (std::size_t ti = 0; ti < src().size(); ++ti)
            if (coarse_box * grow(src()[ti], ghost_nbr))
            {
                host_splits.push_back({ti, dst_box, split_box, coarse_box});
                n_src += src()[ti]().size();
            }
    }
    if (n_src == 0)
        return;

    core::gpu::Vec_t<TileSplit> const splits(host_splits.begin(), host_splits.end());
    core::gpu::Vec_t<Particle_t> staged(n_src * n_kids);
    core::gpu::Vec_t<unsigned long long> counter(1, 0);
    auto const* const sp = splits.data();
    auto* const out      = staged.data();
    auto* const cnt      = counter.data();
    auto view            = *const_cast<Src&>(src);

    Launcher{splits.size()}([=] _PHARE_ALL_FN_() mutable {
        auto const& s = sp[Launcher::block_idx()];
        auto& parts   = view()[s.tile]();
        Particle_t kept[n_kids];
        for (std::size_t i = Launcher::thread_idx(); i < parts.size(); i += Launcher::threads())
        {
            auto const& particle = parts[i];
            if (not isIn(particle, s.coarse_box) or not isIn(fn0(particle), s.split_box))
                continue;
            if (auto const n = fn1(particle, s.dst_box, kept))
            {
                auto const at = Op{*cnt}.add_return_old(n);
                for (std::uint16_t k = 0; k < n; ++k)
                    out[at + k] = kept[k];
            }
        }
    });

    dst.template append_on_device<type>(staged.data(), counter[0]);
#else
    throw std::runtime_error("ParticlesRefiner<AoSTS, GPU_UNIFIED> needs mkn.gpu");
#endif
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSPCTS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                        auto fn1)
{
    ParticlesRefiner<AoSPCTS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}




} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER */
