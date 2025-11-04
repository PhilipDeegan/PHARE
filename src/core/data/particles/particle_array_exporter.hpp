// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER

#include "core/data/particles/exporting/detail/def_exporting.hpp"
#include "core/data/particles/exporting/particles_exporting.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

template<bool in = true, typename Src, typename Box_t>
void delete_particles(Src& src, Box_t const& box)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID     = "delete_particles,";
    [[maybe_unused]] auto constexpr function_id = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

    Exporter{}.template delete_particles<in>(src, box);
}




template<typename Dst, typename Src, typename Box_t>
void move_in_domain(Dst& dst, Src& src, Box_t const& domain_box)
{
    // todo move to dispatcher and optimize
    static_assert(Dst::layout_mode == Src::layout_mode);

    if constexpr (Dst::layout_mode == LayoutMode::AoSTS)
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
    }

    src.sync();
    dst.sync();
}

template<typename Dst, typename Src, typename T, std::size_t dim> //
void move_in_ghost_layer(Dst& dst, Src& src, Box<T, dim> const& domain_box,
                         Box<T, dim> const& ghost_box)
{
    // todo move to dispatcher and optimize
    static_assert(Dst::layout_mode == Src::layout_mode);

    auto const ghost_layer_boxes = ghost_box.remove(domain_box);

    if constexpr (Dst::layout_mode == LayoutMode::AoSTS)
    {
        for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
        {
            auto& src_tile = src()[tidx];
            auto& dst_tile = dst()[tidx];

            auto const tile_ghost_box = grow(src_tile, 1);
            assert(src_tile == dst_tile);

            if (tile_ghost_box * domain_box == tile_ghost_box) // not border tile
                continue;

            auto const in_box = partition_particles(src_tile(), domain_box);
            auto const start  = in_box.size();
            auto const end    = src_tile().size();
            auto const size   = end - start;

            for (std::size_t i = end; i-- > start;)
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
    }

    src.sync();
    dst.sync();
}

template<typename Dst, typename Src, typename T, std::size_t dim, typename Boxes>
void move_in_ghost_layer(Dst& dst, Src& src, Box<T, dim> const& domain_box,
                         Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        move_in_ghost_layer(dst, src, domain_box, gb);
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
