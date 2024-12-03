#ifndef PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP

#include "core/utilities/box/box.hpp"

#include <cassert>
#include <cstdint>
#include <optional>

namespace PHARE::core
{

template<typename TileSet_t, typename Fn>
void traverse_ghost_boundary_tiles(TileSet_t& tileset, Fn fn)
{
    auto const domain_tile_box = shrink(tileset.box(), 2); // bleh

    for (auto& tile : tileset)
        if (auto overlap = domain_tile_box * tile; overlap != tile)
            fn(tile);
}

template<typename TileSet_t, typename Tile_t, typename Fn>
void traverse_tile_neighbours(TileSet_t& tileset, Tile_t& tile, Fn fn)
{
    // 19 = (3**3) - 1 - 7
    auto constexpr backlinks = std::array<std::array<int, 3>, 19>{{
        {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1},
        {-1, 1, 0},   {-1, 1, 1},  {0, -1, -1}, {0, -1, 0},  {0, -1, 1}, {0, 0, -1}, {0, 1, -1},
        {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 1, -1},
    }};

    auto const neighbour_point = [&](auto point) -> std::optional<Point<std::uint32_t, 3>> {
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            if (point[i] < 0)
                point[i] = tile.lower[i] - 1;
            else if (point[i] > 0)
                point[i] = tile.upper[i] + 1;
            else
                point[i] = tile.lower[i];
        }

        if (!isIn(point, tileset.box()))
            return std::nullopt;

        return (point - tileset.box().lower).as_unsigned();
    };

    for (auto const& link : backlinks)
        if (auto const point = neighbour_point(Point{link}))
            fn(*tileset.at(*point));

    for (auto* link : tile.links())
        if (link)
            fn(*link);
}

} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP*/
