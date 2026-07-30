// IWYU pragma: private, include "amr/data/particles/refining/particles_refining.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER


#include "core/def.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_selector.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/refining/detail/def_refining.hpp"

#include <SAMRAI/hier/BoxContainer.h>


namespace PHARE::amr
{

using enum core::LayoutMode;
using enum AllocatorMode;


// for Ghost, particles whose offspring only ever land in the interior dst_amr_box never
// reach the ghost shell, so exclude them via inclusion-exclusion.
template<auto type, typename Src, typename Box_t>
std::size_t relevant_coarse_count(Src const& src, Box_t const& coarseDstBox,
                                  Box_t const& dst_amr_box)
{
    auto const total = core::count_particles(src, coarseDstBox);
    if constexpr (type == core::ParticleType::Ghost)
    {
        auto const coarseAmrBox = coarsen_box(dst_amr_box);
        if (auto const overlap = coarseDstBox * coarseAmrBox)
            return total - core::count_particles(src, *overlap);
    }
    return total;
}

// every one of count's nbRefinedPart children lands somewhere in the box coarseDstBox was
// derived from (that's the whole box, not one fine cell) - no /2**dim here, that would turn
// a box total into a per-cell average.
template<auto options, typename Src>
constexpr std::size_t scale_by_density(std::size_t const count, double const factor)
{
    return static_cast<std::size_t>(count * options.nbRefinedPart * factor);
}

template<auto type, auto options, typename Src, typename Box_t>
std::size_t reserve_flat_count(Src const& src, Box_t const& dst_box, Box_t const& dst_amr_box,
                               double const factor)
{
    // Domain children are clipped to the exact domain box during streaming
    // (ParticlesRefining::forBoxes' per_particle partitions against domainBox). Level-ghost
    // overlaps always sit on top of coarser *domain* cells (AMR nesting guarantees a fine
    // level's ghost region is covered by the coarse level's interior, never its ghost layer),
    // so neither case needs the particle_ghost_width margin here.
    auto const coarseDstBox = coarsen_box(dst_box);
    auto const count        = relevant_coarse_count<type>(src, coarseDstBox, dst_amr_box);
    return scale_by_density<options, Src>(count, factor);
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void reserve_flat(Src const& src, Dst& dst, Box_t const& dst_box, Box_t const& dst_amr_box,
                  double const factor)
{
    dst.reserve(dst.size() + reserve_flat_count<type, options>(src, dst_box, dst_amr_box, factor));
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void reserve_flat(Src const& src, Dst& dst, SAMRAI::hier::BoxContainer const& dst_boxes,
                  Box_t const& dst_amr_box)
{
    std::size_t n = 0;
    for (auto const& samrai_box : dst_boxes)
    {
        n += reserve_flat_count<type, options>(src, phare_box_from<Src::dimension>(samrai_box),
                                               dst_amr_box, .9);

        // grow box by particle_ghost_width
        // remove samrai_box from new bigger box
        // operate across remainder at lower weight
        SAMRAI::hier::BoxContainer ghost_boxes{grow(samrai_box, options.particle_ghost_width)};
        ghost_boxes.removeIntersections(samrai_box);
        for (auto const& border_box : ghost_boxes)
            n += reserve_flat_count<type, options>(src, phare_box_from<Src::dimension>(border_box),
                                                   dst_amr_box, .2);
    }

    dst.reserve(dst.size() + n);
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void stream_split_particles(Src const& src, Dst& dst, Box_t const& dst_box, auto fn0, auto fn1)
{
    std::uint16_t constexpr static N       = 256;
    static constexpr auto base_layout_type = core::base_layout_type<Src>();
    static constexpr auto array_opts
        = Src::options.with_storage(core::StorageMode::ARRAY).with_layout(base_layout_type);
    static constexpr auto array_type_opts
        = core::ParticleArrayTypeOptions<array_opts, base_layout_type,
                                         core::StorageMode::ARRAY>::FROM(Src::options, N);
    using ArrayParticleArray = core::ParticleArrayResolver<array_opts, array_type_opts>::value_type;

    // Resolved the same way as ArrayParticleArray above (not Src::Span_t): that alias
    // only exists on the flat AoS/SoA/SoAVX resolved classes, not e.g. PerCellVector
    // (AoSPC), which is what per_tile_particles resolves to for AoSPCTS.
    static constexpr auto span_opts
        = Src::options.with_storage(core::StorageMode::SPAN).with_layout(base_layout_type);
    static constexpr auto span_type_opts
        = core::ParticleArrayTypeOptions<span_opts, base_layout_type,
                                         core::StorageMode::SPAN>::FROM(Src::options);
    using SpanParticleArray = core::ParticleArrayResolver<span_opts, span_type_opts>::value_type;

    // particle_ghost_width == Splitter::maxCellDistanceFromSplit() (see splitter.hpp) - a
    // coarse particle up to that many fine cells outside dst_box can still have a child land
    // inside it, so this margin is load-bearing for correctness, unlike reserve_flat_count's.
    auto const splitBox     = grow(dst_box, options.particle_ghost_width);
    auto const coarseDstBox = coarsen_box(splitBox);

    std::uint16_t big_buffer_cnt = 0;
    ArrayParticleArray big_buffer;
    SpanParticleArray span{big_buffer};

    auto const send = [&]() {
        assert(big_buffer_cnt <= N);
        span.resize(big_buffer_cnt);
        append_particles<type>(span, dst);
        big_buffer_cnt = 0;
    };

    for (auto const& particle : src)
    {
        if (not isIn(particle, coarseDstBox))
            continue;

        if (not isIn(fn0(particle), splitBox))
            continue;

        auto&& [p_count, buffer] = fn1(particle, dst_box);
        if (big_buffer_cnt + p_count > N)
            send();

        std::copy(buffer.data(), buffer.data() + p_count, big_buffer.data() + big_buffer_cnt);
        big_buffer_cnt += p_count;
    }

    if (big_buffer_cnt)
        send();
}


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
                                                        dst_amr_box);

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
                                                        dst_amr_box);
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



template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSPCTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // per-tile destination is per-cell storage, with no single flat buffer for reserve()
    // to size - see the shrink_to_fit_in comment below.
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;

    auto const& dstBox   = dst.box();
    auto const ghost_nbr = (dst.ghost_box().shape()[0] - dst.box().shape()[0]) / 2;

    using OptBox_t = std::optional<core::Box<int, Src::dimension>>;

    // Per-cell reserve: each dst cell's own coarse-parent cell count (not a box-wide average)
    // - uneven patterns (3D's 6-particle Pink) still spread unevenly across a coarse cell's
    // children, but that's a per-cell rounding slop, not the systematic box-total error a
    // uniform box-average would introduce.
    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const& splitBox    = dst_box;
        auto const coarseDstBox = coarsen_box(dst_box);

        auto const dst_tile_overlap = [&](auto& dst_tile) -> OptBox_t {
            if constexpr (type == core::ParticleType::Domain)
            {
                if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                    return dst_box * dst_tile;
                return std::nullopt;
            }
            else
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
                    return dst_box * excl_box;
                return std::nullopt;
            }
        };

        // reserve pass: one .reserve() call per dst cell, summed across every overlapping
        // src_tile first - the previous version called .reserve() once per (src_tile, cell)
        // pair, multiplying the call count by however many src_tiles touched a cell.
        for (auto& dst_tile : dst())
            if (auto const overlap = dst_tile_overlap(dst_tile))
            {
                auto& dst_pc = dst_tile();
                for (auto const& bix : *overlap)
                {
                    auto const coarseCellPt
                        = coarsen_box(core::Box<int, Src::dimension>{*bix, *bix}).lower;
                    std::size_t coarse_count = 0;
                    for (auto const& src_tile : src())
                    {
                        if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                            continue;
                        auto& src_pc = src_tile();
                        if (isIn(coarseCellPt, src_pc.ghost_box()))
                            coarse_count += src_pc.size(src_pc.local_cell(coarseCellPt));
                    }
                    dst_pc(dst_pc.local_cell(*bix))
                        .reserve(scale_by_density<options, Src>(
                            coarse_count, 1. / (std::size_t{1} << Src::dimension)));
                }
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

    dst.sync();
}



template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSTS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                      auto fn1)
{
    ParticlesRefiner<AoSTS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
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
