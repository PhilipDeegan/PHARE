
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER


#include "core/def.hpp"
#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"


#include <vector>
#include <unordered_map>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>


namespace PHARE::core
{

namespace detail
{
    // tiles of a tile set whose particles can be outside of box: their tile box grown by reach
    //  isn't inside box. Interior tiles have nothing to export
    template<typename Particles, typename Box_t>
    auto tiles_not_inside(Particles const& ps, Box_t const& box, int const reach)
    {
        std::vector<std::size_t> tiles;
        for (std::size_t ti = 0; ti < ps().size(); ++ti)
        {
            auto const& tile = ps()[ti];
            if (tile().size() == 0)
                continue;
            auto const gbox = grow(*tile, reach);
            if (auto const overlap = gbox * box; !overlap or *overlap != gbox)
                tiles.push_back(ti);
        }
        return tiles;
    }

    // tiles of a tile set whose box grown by reach touches any of boxes
    template<typename Particles, typename Boxes>
    auto tiles_touching(Particles const& ps, Boxes const& boxes, int const reach)
    {
        std::vector<std::size_t> tiles;
        for (std::size_t ti = 0; ti < ps().size(); ++ti)
        {
            auto const& tile = ps()[ti];
            if (tile().size() == 0)
                continue;
            auto const gbox = grow(*tile, reach);
            for (auto const& box : boxes)
                if (gbox * box)
                {
                    tiles.push_back(ti);
                    break;
                }
        }
        return tiles;
    }

    // how far from its tile box a tile's particles can be: the ghost layer width, plus the
    //  cell a particle can move past it before the exporters run (move_check leaves particles
    //  leaving the ghost box in place)
    template<typename Particles>
    int tile_reach(Particles const& ps)
    {
        return (ps.ghost_box().shape()[0] - ps.box().shape()[0]) / 2 + 1;
    }

} // namespace detail


using enum LayoutMode;
using enum AllocatorMode;


template<>
template<typename Src, typename Dst, typename Box_t>
void ParticlesExporter<AoSTS, CPU>::move_particles(Src& src, Dst& dst, Box_t const& box,
                                                   std::size_t const growby)
{
    throw std::runtime_error("do not use");

    static_assert(Src::layout_mode == Dst::layout_mode);
    assert(src().size() == dst().size());

    auto const old_size = src.size();
    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile     = src()[tidx];
        auto const src_box = grow(src_tile, growby);
        if (!src_tile().size() or !(box * src_box))
            continue;

        auto& dst_tile        = dst()[tidx];
        auto const not_in_box = partition_particles_not_in(src_tile(), box);
        auto const start      = not_in_box.size();
        auto const end        = src_tile().size();
        auto const size       = end - start;

        PHARE_DEBUG_DO({
            for (std::size_t i = 0; i < start; ++i)
            {
                assert(not isIn(src_tile()[i], box));
                assert(not isIn(src_tile()[i], box));
            }
            for (std::size_t i = start; i < end; ++i)
            {
                assert(isIn(src_tile()[i], box));
            }
            for (std::size_t i = 0; i < dst_tile().size(); ++i)
            {
                assert(isIn(dst_tile()[i], dst_tile));
            }
        })

        if (!size)
            continue;

        append_particles<ParticleType::Domain>(src_tile().view(start, size), dst_tile());

        PHARE_DEBUG_DO({
            for (std::size_t i = 0; i < start; ++i)
            {
                assert(not isIn(src_tile()[i], box));
            }
            for (std::size_t i = start; i < end; ++i)
            {
                assert(isIn(src_tile()[i], box));
            }
            for (std::size_t i = 0; i < dst_tile().size(); ++i)
            {
                assert(isIn(dst_tile()[i], dst_tile));
            }
        })


        if (old_size)
        {
            auto const& p = src_tile()[0];
            assert(start);
        }

        src_tile().resize(start);
    }

    src.template on_appended<ParticleType::Domain>();
    dst.template on_appended<ParticleType::Domain>();

    auto const new_size = src.size();
    if (old_size)
    {
        assert(new_size);
    }
}


template<> // slow
template<typename Src, typename Dst, typename Box_t>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_particles(Src& src, Dst& dst, Box_t const& box,
                                                           std::size_t const growby)
{
    ParticlesExporter<AoSTS, CPU>{}.move_particles(src, dst, box, growby);
}


template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                   Box<int, dim> const& domain_box)
{
    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile = src()[tidx];
        auto& dst_tile = dst()[tidx];

        auto const tile_ghost_box = grow(src_tile, 1);
        assert(src_tile == dst_tile);

        if (tile_ghost_box * domain_box == tile_ghost_box) // not border tile
            continue;

        auto const end = src_tile().size();

        for (std::size_t i = end; i-- > 0;)
        {
            auto const& p = src_tile()[i];
            if (isIn(p, domain_box))
            {
                dst.push_back(p);
                src_tile().vector().erase(src_tile().vector().begin() + i);
            }
        }
    }

    src.template on_appended<ParticleType::Domain>();
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                        Box<int, dim> const& domain_box,
                                                        Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);

    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile = src()[tidx];
        auto& dst_tile = dst()[tidx];

        auto const tile_ghost_box = grow(src_tile, 1);
        assert(src_tile == dst_tile);

        if (tile_ghost_box * domain_box == tile_ghost_box) // not border tile
            continue;

        auto const end = src_tile().size();
        for (std::size_t i = end; i-- > 0;)
        {
            auto const& p = src_tile()[i];
            for (auto const& ghost_layer : ghost_layer_boxes)
            {
                if (isIn(p, ghost_layer))
                {
                    dst_tile().push_back(p);
                    src_tile().vector().erase(src_tile().vector().begin() + i);
                    break;
                }
            }
        }
    }

    src.template on_appended<ParticleType::Domain>();
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                        Box<int, dim> const& domain_box,
                                                        Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


// GPU: one block per tile that can hold such particles, see TileSetVector::move_if_on_device
template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_domain(Dst& dst, Src& src,
                                                           Box<int, dim> const& domain_box)
{
    src.move_if_on_device(dst, detail::tiles_not_inside(src, domain_box, detail::tile_reach(src)),
                          [=] _PHARE_ALL_FN_(auto const& p) { return isIn(p, domain_box); });
}

// GPU: the ghost layer is ghost_box minus domain_box, into the dst tile owning each cell
template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                Box<int, dim> const& domain_box,
                                                                Box<int, dim> const& ghost_box)
{
    src.move_if_on_device(dst, detail::tiles_not_inside(src, domain_box, detail::tile_reach(src)),
                          [=] _PHARE_ALL_FN_(auto const& p) {
                              return isIn(p, ghost_box) and !isIn(p, domain_box);
                          });
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                Box<int, dim> const& domain_box,
                                                                Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return isIn(p.iCell(), box); });
    src.erase(kept, src.end());
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return isIn(cell, boxes); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return isIn(cell, box); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return isIn(p.iCell(), boxes); });
    src.erase(kept, src.end());
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return !isIn(p.iCell(), boxes); });
    src.erase(kept, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                 Box<int, dim> const& domain_box)
{
    auto const moved = std::partition(src.begin(), src.end(),
                                      [&](auto const& p) { return !isIn(p.iCell(), domain_box); });
    std::copy(moved, src.end(), std::back_inserter(dst));
    src.erase(moved, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                      Box<int, dim> const& domain_box,
                                                      Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);
    auto const moved             = std::partition(
        src.begin(), src.end(), [&](auto const& p) { return !isIn(p.iCell(), ghost_layer_boxes); });
    std::copy(moved, src.end(), std::back_inserter(dst));
    src.erase(moved, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                      Box<int, dim> const& domain_box,
                                                      Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return !isIn(cell, boxes); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, GPU_UNIFIED>::delete_particles_not_in(Src& src,
                                                                        Boxes const& boxes)
{
    throw std::runtime_error("todo");
}

template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto const ghost_nbr = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        if (auto const tile_gbox = grow(*tile, ghost_nbr); tile_gbox * box)
        {
            auto const range = partition_particles(tile(), box);
            auto const resiz = range.size();

            assert(resiz <= tile().size());
            tile().resize(resiz);
        }
    }

    src.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto const ghost_nbr = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        if (auto const tile_gbox = grow(*tile, ghost_nbr); tile_gbox * src.box() != tile_gbox)
        {                                                           // only patch border tiles
            auto const ranges = partition_particles(tile(), boxes); // what
            auto const resiz  = sum_from(ranges, [](auto const& r) { return r.size(); });

            assert(resiz <= tile().size());
            tile().resize(resiz);
        }
    }

    src.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto const ghost_nbr = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& tile : src())
        for (auto const& box : boxes)
            if (tile().size() and grow(*tile, ghost_nbr) * box)
                tile().resize(partition_particles_not_in(tile(), box).size());

    src.template on_appended<ParticleType::Domain>();
}



// AoSCMTS: one flat cell-mapped array per tile, see MappedTileSetVector::erase_if
template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    src.erase_if([&](auto const& p) { return !isIn(p.iCell(), box); });
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    src.erase_if([&](auto const& p) { return !isIn(p.iCell(), boxes); });
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), boxes); });
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                     Box<int, dim> const& domain_box)
{
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), domain_box); },
                 [&](auto const& p) { dst.push_back(p); });
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), ghost_layer_boxes); },
                 [&](auto const& p) { dst.push_back(p); });
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


// GPU: one block per tile reaching outside box, see TileSetVector::erase_if_on_device
template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles_not_in(Src& src,
                                                                    Box<int, dim> const& box)
{
    src.erase_if_on_device(detail::tiles_not_inside(src, box, detail::tile_reach(src)),
                           [=] _PHARE_ALL_FN_(auto const& p) { return !isIn(p, box); });
}

// boxes go to managed memory for the device predicate
template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
#if PHARE_HAVE_MKN_GPU
    using Box_t = std::decay_t<decltype(*std::begin(boxes))>;
    gpu::Vec_t<Box_t> const dboxes(std::begin(boxes), std::end(boxes));
    auto const* const bp = dboxes.data();
    auto const nb        = dboxes.size();
    src.erase_if_on_device(detail::tiles_not_inside(src, src.box(), detail::tile_reach(src)),
                           [=] _PHARE_ALL_FN_(auto const& p) {
                               for (std::size_t i = 0; i < nb; ++i)
                                   if (isIn(p, bp[i]))
                                       return false;
                               return true;
                           });
#else
    ParticlesExporter<AoSTS, CPU>{}.delete_particles_not_in(src, boxes);
#endif
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles_in(Src& src, Boxes const& boxes)
{
#if PHARE_HAVE_MKN_GPU
    using Box_t = std::decay_t<decltype(*std::begin(boxes))>;
    gpu::Vec_t<Box_t> const dboxes(std::begin(boxes), std::end(boxes));
    auto const* const bp = dboxes.data();
    auto const nb        = dboxes.size();
    src.erase_if_on_device(detail::tiles_touching(src, boxes, detail::tile_reach(src)),
                           [=] _PHARE_ALL_FN_(auto const& p) {
                               for (std::size_t i = 0; i < nb; ++i)
                                   if (isIn(p, bp[i]))
                                       return true;
                               return false;
                           });
#else
    ParticlesExporter<AoSTS, CPU>{}.delete_particles_in(src, boxes);
#endif
}


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for_each_point(cps.local_box(gb), [&](auto const& bix) {
            auto& cell_parts = cps(bix);
            if (cell_parts.size() == 0)
                return;
            std::array<int, dim> gcell;
            for (std::size_t d = 0; d < dim; ++d)
                gcell[d] = gb.lower[d] + static_cast<int>(bix[d]);
            if (!isIn(gcell, box))
                cell_parts.clear();
        });
        cps.template on_appended<ParticleType::Ghost>();
    }
    src.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSPCTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto constexpr dim = Src::dimension;
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for_each_point(cps.local_box(gb), [&](auto const& bix) {
            auto& cell_parts = cps(bix);
            if (cell_parts.size() == 0)
                return;
            std::array<int, dim> gcell;
            for (std::size_t d = 0; d < dim; ++d)
                gcell[d] = gb.lower[d] + static_cast<int>(bix[d]);
            bool in_any = false;
            for (auto const& b : boxes)
                if (isIn(gcell, b))
                {
                    in_any = true;
                    break;
                }
            if (!in_any)
                cell_parts.clear();
        });
        cps.template on_appended<ParticleType::Ghost>();
    }
    src.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSPCTS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto constexpr dim = Src::dimension;
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for_each_point(cps.local_box(gb), [&](auto const& bix) {
            auto& cell_parts = cps(bix);
            if (cell_parts.size() == 0)
                return;
            std::array<int, dim> gcell;
            for (std::size_t d = 0; d < dim; ++d)
                gcell[d] = gb.lower[d] + static_cast<int>(bix[d]);
            if (isIn(gcell, boxes))
                cell_parts.clear();
        });
        cps.template on_appended<ParticleType::Ghost>();
    }
    src.template on_appended<ParticleType::Domain>();
}


// AoSPCTS domain <-> patch-ghost/level-ghost particle exchange: per-cell aware, mirrors
// delete_particles_not_in<AoSPCTS, CPU> traversal. `dst.push_back()` resolves the
// receiving tile+cell from the particle's own iCell(), so src/dst tiles need not align.
template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                     Box<int, dim> const& domain_box)
{
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for_each_point(cps.local_box(gb), [&](auto const& bix) {
            auto& cell_parts = cps(bix);
            for (std::size_t i = cell_parts.size(); i-- > 0;)
            {
                if (isIn(cell_parts.iCell(i), domain_box))
                {
                    dst.push_back(cell_parts[i]);
                    cell_parts.assign(cell_parts.size() - 1, i);
                    cell_parts.pop_back();
                }
            }
        });
        cps.template on_appended<ParticleType::Ghost>();
    }
    src.template on_appended<ParticleType::Domain>();
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);

    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for_each_point(cps.local_box(gb), [&](auto const& bix) {
            auto& cell_parts = cps(bix);
            for (std::size_t i = cell_parts.size(); i-- > 0;)
            {
                auto const& icell = cell_parts.iCell(i);
                for (auto const& ghost_layer : ghost_layer_boxes)
                {
                    if (isIn(icell, ghost_layer))
                    {
                        dst.push_back(cell_parts[i]);
                        cell_parts.assign(cell_parts.size() - 1, i);
                        cell_parts.pop_back();
                        break;
                    }
                }
            }
        });
        cps.template on_appended<ParticleType::Ghost>();
    }
    src.template on_appended<ParticleType::Domain>();
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSPCTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


// AoSPCTS GPU_UNIFIED: a cell's bucket only holds that cell's particles, so the exporters
//  decide per cell. Deleting clears buckets (no particle data is touched): the CPU versions.
//  Moving copies whole buckets on the device, see detail::pc_gpu::move_cells

namespace detail::pc_gpu
{
    template<typename Particle_t>
    struct CellCopy
    {
        Particle_t const* src;
        Particle_t* dst;
        std::size_t n;
    };

    // moves the buckets of the src cells that pred(amr cell) selects into the dst bucket of the
    //  same cell (in the dst tile owning it): buckets are sized on the host, particles copied on
    //  the device, one block per cell
    template<typename Dst, typename Src>
    void move_cells(Dst& dst, Src& src, auto&& pred)
    {
#if PHARE_HAVE_MKN_GPU
        auto constexpr dim = Src::dimension;
        using Particle_t   = typename Src::Particle_t;
        using Launcher     = gpu::BlockPerLauncher<true>;

        struct Job
        {
            void* src_tile;
            void* dst_tile;
            void* src_bucket;
            void* dst_bucket;
            std::size_t n;
        };
        std::vector<Job> jobs;
        std::unordered_map<void*, std::size_t> added; // per dst bucket

        for (auto& tile : src())
        {
            auto& cps      = tile();
            auto const& gb = cps.ghost_box();
            for_each_point(cps.local_box(gb), [&](auto const& bix) {
                auto& bucket = cps(bix);
                if (bucket.size() == 0)
                    return;
                std::array<int, dim> amr;
                for (std::size_t d = 0; d < dim; ++d)
                    amr[d] = gb.lower[d] + static_cast<int>(bix[d]);
                if (!pred(amr))
                    return;
                auto* dst_tile = dst().at(Point<int, dim>{amr});
                assert(dst_tile);
                auto& dcps = (*dst_tile)();
                auto& dst_bucket = dcps(dcps.local_cell(amr));
                jobs.push_back({&cps, &dcps, &bucket, &dst_bucket, bucket.size()});
                added[&dst_bucket] += bucket.size();
            });
        }
        if (jobs.empty())
            return;

        using PC_t     = std::decay_t<decltype(src()[0]())>;
        using Bucket_t = std::decay_t<decltype(src()[0]()(std::array<std::uint32_t, dim>{}))>;

        std::unordered_map<void*, std::size_t> cursor; // dst buckets grow once, before copying
        for (auto const& [bucket, n] : added)
        {
            auto& b         = *static_cast<Bucket_t*>(bucket);
            cursor[bucket]  = b.size();
            PC_t::resize(b, b.size() + n);
        }

        std::vector<CellCopy<Particle_t>> host_copies;
        host_copies.reserve(jobs.size());
        for (auto const& job : jobs)
        {
            auto& sb   = *static_cast<Bucket_t*>(job.src_bucket);
            auto& db   = *static_cast<Bucket_t*>(job.dst_bucket);
            auto& at   = cursor[job.dst_bucket];
            host_copies.push_back({sb.data(), db.data() + at, job.n});
            at += job.n;
        }

        gpu::Vec_t<CellCopy<Particle_t>> const copies(host_copies.begin(), host_copies.end());
        auto const* const cp = copies.data();
        Launcher{copies.size()}([=] _PHARE_ALL_FN_() {
            auto const& c = cp[Launcher::block_idx()];
            for (std::size_t i = Launcher::thread_idx(); i < c.n; i += Launcher::threads())
                c.dst[i] = c.src[i];
        });

        std::unordered_map<void*, bool> tiles; // src and dst tiles to refresh
        for (auto const& job : jobs)
        {
            static_cast<Bucket_t*>(job.src_bucket)->clear();
            tiles[job.src_tile] = true;
            tiles[job.dst_tile] = true;
        }
        for (auto const& [tile, _] : tiles)
            static_cast<PC_t*>(tile)->template on_appended<ParticleType::Ghost>();
        src.template on_appended<ParticleType::Domain>();
        dst.template on_appended<ParticleType::Domain>();
#else
        throw std::runtime_error("AoSPCTS GPU particle moves need mkn.gpu");
#endif
    }

} // namespace detail::pc_gpu


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::delete_particles_not_in(Src& src,
                                                                      Box<int, dim> const& box)
{
    ParticlesExporter<AoSPCTS, CPU>{}.delete_particles_not_in(src, box);
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    ParticlesExporter<AoSPCTS, CPU>{}.delete_particles_not_in(src, boxes);
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::delete_particles_in(Src& src, Boxes const& boxes)
{
    ParticlesExporter<AoSPCTS, CPU>{}.delete_particles_in(src, boxes);
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::move_in_domain(Dst& dst, Src& src,
                                                             Box<int, dim> const& domain_box)
{
    detail::pc_gpu::move_cells(dst, src,
                               [&](auto const& amr) { return isIn(amr, domain_box); });
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                  Box<int, dim> const& domain_box,
                                                                  Box<int, dim> const& ghost_box)
{
    detail::pc_gpu::move_cells(dst, src, [&](auto const& amr) {
        return isIn(amr, ghost_box) and !isIn(amr, domain_box);
    });
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSPCTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                  Box<int, dim> const& domain_box,
                                                                  Boxes const& ghost_boxes)
{
    detail::pc_gpu::move_cells(dst, src, [&](auto const& amr) {
        return isIn(amr, ghost_boxes) and !isIn(amr, domain_box);
    });
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER */
