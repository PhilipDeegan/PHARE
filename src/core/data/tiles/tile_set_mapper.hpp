#ifndef PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP

#include "core/data/tiles/tile_set.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include <cassert>
#include <cstdint>

namespace PHARE::core
{

template<typename TileSet_t>
struct ATiler
{
    using Tile_t                    = TileSet_t::value_type;
    using Box_t                     = TileSet_t::Box_t;
    static auto constexpr dimension = TileSet_t::dimension;

    TileSet_t& tile_set;
    Box_t const& box                                    = tile_set.box();
    std::array<std::uint32_t, dimension> const& shape   = tile_set.shape();
    std::array<std::size_t, dimension> const& tile_size = tile_set.tile_size();
};

template<typename TileSet_t, int impl = 0>
struct Tiler;


template<typename TileSet_t>
struct Tiler<TileSet_t, 0> : public ATiler<TileSet_t>
{
    template<typename... Args>
    void f(Args&&... args);
};

template<typename TileSet_t>
template<typename... Args>
void Tiler<TileSet_t, 0>::f(Args&&... args)
{
    using Tile_t                    = TileSet_t::value_type;
    using Box_t                     = TileSet_t::Box_t;
    static auto constexpr dimension = TileSet_t::dimension;

    auto const& box       = this->box;
    auto const& shape     = this->shape;
    auto const& tile_size = this->tile_size;

    auto const size_me = [&](auto dim, auto idx) {
        if (idx == shape[dim] - 1)
        {
            auto const remain = box.shape()[dim] % tile_size[dim];
            return (remain == 0) ? tile_size[dim] : remain;
        }
        else
            return tile_size[dim];
    };

    auto& tiles = this->tile_set();

    for (auto ix = 0u; ix < shape[0]; ++ix)
    {
        Box_t tile;
        if constexpr (dimension == 1)
        {
            // -1 because upper is included
            tile.lower[0] = box.lower[0] + ix * tile_size[0];
            tile.upper[0] = tile.lower[0] + size_me(0, ix) - 1;
            tiles.emplace_back(tile, args...);
        }
        else
        {
            for (auto iy = 0u; iy < shape[1]; ++iy)
            {
                if constexpr (dimension == 2)
                {
                    // auto const i  = ix * shape[1] + iy;
                    tile.lower[0] = box.lower[0] + ix * tile_size[0];
                    tile.upper[0] = tile.lower[0] + size_me(0, ix) - 1;
                    tile.lower[1] = box.lower[1] + iy * tile_size[1];
                    tile.upper[1] = tile.lower[1] + size_me(1, iy) - 1;
                    tiles.emplace_back(tile, args...);
                }
                else
                {
                    for (auto iz = 0u; iz < shape[2]; ++iz)
                    {
                        // auto const i  = ix * shape[1] * shape[2] + shape[2] * iy + iz;
                        tile.lower[0] = box.lower[0] + ix * tile_size[0];
                        tile.upper[0] = tile.lower[0] + size_me(0, ix) - 1;
                        tile.lower[1] = box.lower[1] + iy * tile_size[1];
                        tile.upper[1] = tile.lower[1] + size_me(1, iy) - 1;
                        tile.lower[2] = box.lower[2] + iz * tile_size[2];
                        tile.upper[2] = tile.lower[2] + size_me(2, iz) - 1;
                        tiles.emplace_back(tile, args...);
                    }
                }
            }
        }
    }
}



template<typename TileSet_t>
struct Embigginer
{
    auto constexpr static dim = TileSet_t::dimension;

    Embigginer(ATiler<TileSet_t> const& tiler_)
        : tiler{tiler_}
    {
        //     self.count = n_tiles if rem==n_tiles  else rem if int(n_tiles /2)>rem else 1
        //     self.add = 1 if self.count == rem else int(rem / self.count)

        //     self.startx = int(n_tiles  / 2) if self.count == 1 else int((n_tiles -1)  / 2)
        //     print(self.count, self.add, self.startx)
    }

    void doDim(int const d, int const xyz)
    {
        add = tiler.tile_size[d] - 1;
        if (xyz >= start[d] and nd[d] < count[d])
        {
            dd[d] += add;
            nd[d] += 1;
        }
    }

    ATiler<TileSet_t> const& tiler;
    std::array<int, dim> count         = ConstArray<int, dim>();
    std::array<int, dim> add           = ConstArray<int, dim>();
    std::array<int, dim> dd            = ConstArray<int, dim>();
    std::array<int, dim> nd            = ConstArray<int, dim>();
    std::array<int, dim> start         = ConstArray<int, dim>();
    std::array<int, dim> tiles_per_dim = ConstArray<int, dim>();
};

template<typename TileSet_t>
struct Tiler<TileSet_t, 1> : public ATiler<TileSet_t>
{
    template<typename... Args>
    void f(Args&&... args);
};

template<typename TileSet_t>
template<typename... Args>
void Tiler<TileSet_t, 1>::f(Args&&... args)
{
    using Embigginer_t        = Embigginer<TileSet_t>;
    auto constexpr static dim = TileSet_t::dimension;
    Embigginer_t e{*this};

    // assuming 3d for now
    if constexpr (dim == 1)
    {
        throw "fail";
    }
    if constexpr (dim == 2)
    {
        throw "fail";
    }
    else
    {
        std::array<int, dim> add  = ConstArray<int, dim>();
        std::array<int, dim> cell = ConstArray<int, dim>();

        for (std::uint8_t i = 0; i < e.tiles_per_dim[0]; ++i)
        {
            e.doDim(0, i);
            for (std::uint8_t j = 0; j < e.tiles_per_dim[1]; ++j)
            {
                e.doDim(1, j);

                for (std::uint8_t k = 0; k < e.tiles_per_dim[2]; ++k)
                {
                    e.doDim(2, k);

                    int xlo = cell[0], xhi = cell[0] + e.dd[0] + this->tile_size[0] - 1;
                    int ylo = cell[1], yhi = cell[1] + e.dd[1] + this->tile_size[1] - 1;
                    int zlo = cell[2], zhi = cell[2] + e.dd[2] + this->tile_size[2] - 1;

                    Point<int, dim> const lo{xlo, ylo, zlo};
                    Point<int, dim> const hi{xhi, yhi, zhi};

                    //
                    cell[2] += this->tile_size[2] + e.dd[2];
                    e.z = 0;
                }
                cell[1] += this->tile_size[1] + e.dd[1];
                cell[2] = 0;
                e.dd[1] = 0;
                e.nd[2] = 0;
            }

            cell[0] += this->tile_size[0] + e.x;
            cell[1] = 0;
            cell[2] = 0;
            e.dd[0] = 0;
            e.nd[1] = 0;
        }
    }
}

template<std::uint8_t impl = 0, typename TileSet_t, typename... Args>
void tile_set_make_tiles(TileSet_t& tile_set, Args&&... args)
{
    Tiler<TileSet_t, impl>{{tile_set}}.f(args...);
}



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP*/
