
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER


#include "core/def.hpp"
#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"


#include <stdexcept>


namespace PHARE::core
{

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

    src.sync();
    dst.sync();

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
template<bool in, typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::delete_particles(Src& src, Box<int, dim> const& box)
{
    throw std::runtime_error("todo");
}
template<>
template<bool in, typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles(Src& src, Boxes const& boxes)
{
    throw std::runtime_error("todo");
}

template<>
template<bool in, typename Src, std::size_t dim>
void ParticlesExporter<AoSMapped, CPU>::delete_particles(Src& src, Box<int, dim> const& box)
{
    throw std::runtime_error("todo");
}
template<>
template<bool in, typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, GPU_UNIFIED>::delete_particles(Src& src, Boxes const& boxes)
{
    throw std::runtime_error("todo");
}

template<>
template<bool in, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::delete_particles(Src& src, Box<int, dim> const& box)
{
    static_assert(in == false); // for now

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

    src.sync();
}

template<>
template<bool in, typename Src, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::delete_particles(Src& src, Boxes const& boxes)
{
    static_assert(in == false); // for now

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

    src.sync();
}




template<>
template<bool in, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles(Src& src, Box<int, dim> const& box)
{
    ParticlesExporter<AoSTS, CPU>{}.delete_particles<in>(src, box);
}

template<>
template<bool in, typename Src, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles(Src& src, Boxes const& boxes)
{
    ParticlesExporter<AoSTS, CPU>{}.delete_particles<in>(src, boxes);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER */
