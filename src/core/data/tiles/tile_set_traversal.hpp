#ifndef PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP

#include <cassert>
#include <functional>

namespace PHARE::core
{

template<typename TileSet_t, typename Box_t, typename Fn>
void traverse_tiles(TileSet_t& tileset, Box_t const& box, Fn fn)
{
    using Tile_t = typename TileSet_t::value_type;

    std::function<void(Tile_t&)> const doX = [&](auto& tile) {
        fn(tile);

        if (auto nextX = tile.link(0); nextX and (box * (**nextX)))
            doX(*nextX);
    };
    std::function<void(Tile_t&)> const doY = [&](auto& tile) {
        doX(tile);

        if constexpr (TileSet_t::dimension > 1)
            if (auto nextY = tile.link(1); nextY and (box * (**nextY)))
                doY(*nextY);
    };
    std::function<void(Tile_t&)> const doZ = [&](auto& tile) {
        doY(tile);

        if constexpr (TileSet_t::dimension == 3)
            if (auto nextZ = tile.link(3); nextZ and (box * (**nextZ)))
                doZ(*nextZ);
    };

    doZ(*tileset.at(box.lower));
}

} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP*/
