#ifndef PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP

#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

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
    Box_t const& box = tile_set.box();
};

template<typename TileSet_t, int impl = 0>
struct Tiler;


template<typename TileSet_t>
struct Tiler<TileSet_t, 0> : public ATiler<TileSet_t>
{
    template<typename... Args>
    auto f(Args&&... args);

    std::array<std::uint32_t, ATiler<TileSet_t>::dimension> const& shape;
    std::array<std::size_t, ATiler<TileSet_t>::dimension> const& tile_size;
};

template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t, 0>::f(Args&&... args)
{
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
struct EmbigginerBase
{
    auto constexpr static dim = TileSet_t::dimension;

    EmbigginerBase(ATiler<TileSet_t> const& tiler_)
        : tiler{tiler_}
    {
    }

    void reset()
    {
        auto const shape = tiler.box.shape();
        for (std::uint8_t i = 0; i < dim; ++i)
        {
            count[i] = rem_per_dim[i] == tiles_per_dim[i]           ? tiles_per_dim[i]
                       : int(tiles_per_dim[i] / 2) > rem_per_dim[i] ? rem_per_dim[i]
                                                                    : 1;
            add[i]   = count[i] == rem_per_dim[i] ? 1 : int(rem_per_dim[i] / count[i]);
            start[i] = count[i] == tiles_per_dim[i] ? 0
                       : count[i] == 1              ? int(tiles_per_dim[i] / 2)
                                                    : int((tiles_per_dim[i] - 1) / 2);
        }
    }

    void doDim(int const d, int const xyz)
    {
        if (xyz >= start[d] and nd[d] < count[d])
        {
            dd[d] += add[d];
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
    std::array<int, dim> rem_per_dim   = ConstArray<int, dim>();
};


template<typename TileSet_t>
struct Embigginer : public EmbigginerBase<TileSet_t>
{
    auto constexpr static dim = TileSet_t::dimension;

    Embigginer(ATiler<TileSet_t> const& tiler_)
        : EmbigginerBase<TileSet_t>{tiler_}
    {
        auto const shape = this->tiler.box.shape();
        for (std::uint8_t i = 0; i < dim; ++i)
        {
            this->tiles_per_dim[i] = shape[i] / this->tiler.tile_size[i];
            this->rem_per_dim[i]   = shape[i] % this->tiler.tile_size[i];
        }
        this->reset();
    }
};

template<typename TileSet_t>
struct Tiler<TileSet_t, 1> : public ATiler<TileSet_t>
{
    template<typename... Args>
    auto f(Args&&... args);

    std::array<std::uint32_t, ATiler<TileSet_t>::dimension> const& shape = this->tile_set.shape();
    std::array<std::size_t, ATiler<TileSet_t>::dimension> const& tile_size
        = this->tile_set.tile_size();
};

template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t, 1>::f(Args&&... args)
{
    auto constexpr static dim = TileSet_t::dimension;

    using Box_t        = TileSet_t::Box_t;
    using Embigginer_t = Embigginer<TileSet_t>;

    auto& tiles = this->tile_set();

    Embigginer_t e{*this};

    // assuming 3d for now
    if constexpr (dim == 1)
    {
        auto cell = ConstArray<int, dim>();

        for (std::uint8_t i = 0; i < e.tiles_per_dim[0]; ++i)
        {
            e.doDim(0, i);

            int const xlo = cell[0], xhi = cell[0] + e.dd[0] + this->tile_size[0] - 1;

            Point<int, dim> const lo{xlo};
            Point<int, dim> const hi{xhi};
            tiles.emplace_back(Box_t{lo, hi}, args...);


            cell[0] += this->tile_size[0] + e.dd[0];
            e.dd[0] = 0;
        }
    }

    if constexpr (dim == 2)
    {
        auto cell = ConstArray<int, dim>();

        for (std::uint8_t i = 0; i < e.tiles_per_dim[0]; ++i)
        {
            e.doDim(0, i);
            for (std::uint8_t j = 0; j < e.tiles_per_dim[1]; ++j)
            {
                e.doDim(1, j);


                int const xlo = cell[0], xhi = cell[0] + e.dd[0] + this->tile_size[0] - 1;
                int const ylo = cell[1], yhi = cell[1] + e.dd[1] + this->tile_size[1] - 1;

                Point<int, dim> const lo{xlo, ylo};
                Point<int, dim> const hi{xhi, yhi};
                tiles.emplace_back(Box_t{lo, hi}, args...);

                cell[1] += this->tile_size[1] + e.dd[1];
                e.dd[1] = 0;
            }

            cell[0] += this->tile_size[0] + e.dd[0];
            cell[1] = 0;
            e.dd[0] = 0;
            e.nd[1] = 0;
        }
    }
    if constexpr (dim == 3)
    {
        auto cell = ConstArray<int, dim>();

        for (std::uint8_t i = 0; i < e.tiles_per_dim[0]; ++i)
        {
            e.doDim(0, i);
            for (std::uint8_t j = 0; j < e.tiles_per_dim[1]; ++j)
            {
                e.doDim(1, j);

                for (std::uint8_t k = 0; k < e.tiles_per_dim[2]; ++k)
                {
                    e.doDim(2, k);

                    int const xlo = cell[0], xhi = cell[0] + e.dd[0] + this->tile_size[0] - 1;
                    int const ylo = cell[1], yhi = cell[1] + e.dd[1] + this->tile_size[1] - 1;
                    int const zlo = cell[2], zhi = cell[2] + e.dd[2] + this->tile_size[2] - 1;

                    Point<int, dim> const lo{xlo, ylo, zlo};
                    Point<int, dim> const hi{xhi, yhi, zhi};
                    tiles.emplace_back(Box_t{lo, hi}, args...);

                    cell[2] += this->tile_size[2] + e.dd[2];
                    e.dd[2] = 0;
                }
                cell[1] += this->tile_size[1] + e.dd[1];
                cell[2] = 0;
                e.dd[1] = 0;
                e.nd[2] = 0;
            }

            cell[0] += this->tile_size[0] + e.dd[0];
            cell[1] = 0;
            cell[2] = 0;
            e.dd[0] = 0;
            e.nd[1] = 0;
        }
    }
}



template<typename TileSet_t>
struct Tiler<TileSet_t, 2> : public ATiler<TileSet_t>
{
    template<typename... Args>
    auto f(Args&&... args);
};

template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t, 2>::f(Args&&... args)
{
    auto constexpr static dim = TileSet_t::dimension;

    using Box_t       = TileSet_t::Box_t;
    auto const& shape = this->box.shape();
    auto& tiles       = this->tile_set();

    core::Point<int, dim> rem_per_dim{ConstArray<int, dim>()};
    core::Point<int, dim> tiles_per_dim{ConstArray<int, dim>()};

    for_N<dim>([&](auto i) {
        if (shape[i] < 10)
            tiles_per_dim[i] = 1;
        else if (shape[i] == 10)
            tiles_per_dim[i] = 2;
        else
            tiles_per_dim[i] = 3;

        if (tiles_per_dim[i] > 2)
            rem_per_dim[i] = shape[i] - 8;
    });

    std::array<std::array<int, 5>, dim> starts;

    auto setDim = [&](auto d, auto... args) {
        std::uint8_t i = 0;
        ((starts[d][i++] = args), ...);
    };

    auto const doDim = [&](int const d) {
        if (tiles_per_dim[d] == 1)
            setDim(d, shape[d]);
        else if (tiles_per_dim[d] == 2)
            setDim(d, shape[d] / 2, shape[d] / 2);
        else if (tiles_per_dim[d] == 3)
            setDim(d, 4, shape[d] - 8, 4);
        // else if (tiles_per_dim[d] == 4)
        //     setDim(d, 4, (shape[d] - 8) / 2, (shape[d] - 8) / 2, 4);
        // else if (tiles_per_dim[d] == 5)
        //     setDim(d, 4, (shape[d] - 8) / 3, ((shape[d] - 8) / 3) + rem_per_dim[d],
        //            (shape[d] - 8) / 3, 4);
    };

    for (std::uint8_t d = 0; d < dim; ++d)
        doDim(d);


    // PHARE_LOG_LINE_SS(tiles_per_dim);

    auto off = ConstArray<int, dim>();

    using Point_t = Point<int, dim>;

    if constexpr (dim == 1)
        for (std::uint8_t i = 0; i < tiles_per_dim[0]; ++i)
        {
            int const xlo = off[0], xhi = off[0] + starts[0][i] - 1;
            Point_t const lo{xlo};
            Point_t const hi{xhi};
            tiles.emplace_back(Box_t{lo, hi}, args...);
            off[0] = xhi + 1;
        }


    if constexpr (dim == 3)
    {
        using TileFn = std::function<void(Point_t, Point_t, Point_t)>;

        TileFn doX = [&](auto from, auto to, auto by) {
            // PHARE_LOG_LINE_SS(from << " " << to << " " << by);
            while (from < to)
            {
                auto up = from;
                up += by - 1;
                // PHARE_LOG_LINE_SS(from << " " << up << " " << by);

                tiles.emplace_back(Box_t{from, up}, args...);
                from[2] += by[2];
            }
        };
        TileFn doY = [&](auto from, auto to, auto by) {
            while (from < to)
            {
                doX(from, to, by);
                from[1] += by[1];
            }
        };
        TileFn doZ = [&](auto from, auto to, auto by) {
            // PHARE_LOG_LINE_SS(from << " " << to << " " << by);
            while (from < to)
                doY(from, to, by), from[0] += by[0];
        };

        std::array<Box_t, 5> boxes;

        using Embigginer_t = EmbigginerBase<TileSet_t>;

        Embigginer_t e{*this};

        auto cell = ConstArray<int, dim>();

        auto const max_dim = *std::max_element(tiles_per_dim.begin(), tiles_per_dim.end());
        auto const min_dim = *std::min_element(tiles_per_dim.begin(), tiles_per_dim.end());

        // PHARE_LOG_LINE_SS(max_dim);

        if (max_dim < 3 and min_dim < 3)
        {
            for (std::uint8_t i = 0; i < tiles_per_dim[0]; ++i)
            {
                for (std::uint8_t j = 0; j < tiles_per_dim[1]; ++j)
                {
                    for (std::uint8_t k = 0; k < tiles_per_dim[2]; ++k)
                    {
                        int const xlo = cell[0], xhi = cell[0] + e.dd[0] + starts[0][i] - 1;
                        int const ylo = cell[1], yhi = cell[1] + e.dd[1] + starts[1][j] - 1;
                        int const zlo = cell[2], zhi = cell[2] + e.dd[2] + starts[2][k] - 1;

                        Point<int, dim> const lo{xlo, ylo, zlo};
                        Point<int, dim> const hi{xhi, yhi, zhi};
                        tiles.emplace_back(Box_t{lo, hi}, args...);

                        cell[2] = starts[2][k] + e.dd[2];
                        e.dd[2] = 0;
                    }
                    cell[1] = starts[1][j] + e.dd[1];
                    cell[2] = 0;
                    e.dd[1] = 0;
                    e.nd[2] = 0;
                }

                cell[0] = starts[0][i] + e.dd[0];
                cell[1] = 0;
                cell[2] = 0;
                e.dd[0] = 0;
                e.nd[1] = 0;
            }
        }
        else if (min_dim > 2)
        {
            auto blocking = [&](auto box, auto by, auto rm) {
                for (auto const rem : box.remove(rm))
                    doZ(rem.lower, rem.upper, by);
            };


            if (for_N_all<dim>([&](auto i) { return rem_per_dim[i] % 4 == 0; }))
            {
                PHARE_LOG_LINE_SS(rem_per_dim);
                Point<int, dim> const by{starts[0][0], starts[1][0], starts[2][0]};
                Box_t rem{{starts[0][0], starts[1][0], starts[2][0]},
                          {starts[0][0] + starts[0][1] - 1, starts[1][0] + starts[1][1] - 1,
                           starts[2][0] + starts[2][1] - 1}};
                blocking(this->box, by, rem);
                doZ(rem.lower, rem.upper, rem.shape());
            }
            else
            {
                PHARE_LOG_LINE_SS(rem_per_dim);
                Point<int, dim> const mid{for_N_make_array<3>([&](auto i) {
                    // if (rem_per_dim[i] > 7)
                    // {
                    //     auto const mod = rem_per_dim[i] % 4;
                    //     auto const div = rem_per_dim[i] / 4;
                    //     auto const sid = div / 2;
                    //     return sid;
                    //     // PHARE_LOG_LINE_SS(rem_per_dim);
                    //     // return rem_per_dim[i] % 2 == 0 ? 8 : 7;
                    // }
                    return rem_per_dim[i];
                })};

                Point<int, dim> const corner_size{
                    for_N_make_array<3>([&](auto i) { return (shape[i] - mid[i]) / 2; })};

                auto shifts = PHARE::core::for_N<8, PHARE::core::for_N_R_mode::make_array>(
                    [&](auto i) { return PHARE::core::Point<int, dim>{0, 0, 0}; });
                PHARE::core::for_N<dim>([&](auto i) { shifts[i][i] = corner_size[i] + mid[i]; });

                shifts[3] = {shifts[0][0], shifts[1][1], 0};
                shifts[4] = {0, shifts[1][1], shifts[2][2]};
                shifts[5] = {shifts[0][0], 0, shifts[2][2]};
                shifts[6] = {shifts[0][0], shifts[1][1], shifts[2][2]};
                shifts[7] = {0, 0, 0};

                Point<int, dim> const by{4, 4, 4};
                auto const& box = this->box;
                Box_t corner{box.lower, box.lower + corner_size - 1};

                for (auto const shifter : shifts)
                {
                    auto const corn = shift(corner, shifter);
                    doZ(corn.lower, corn.upper, by);
                }

                Box_t const middle{corner.upper + 1, corner.upper + mid};

                PHARE_LOG_LINE_SS(middle);
                PHARE_LOG_LINE_SS(corner);
                PHARE_LOG_LINE_SS(shift(corner, shifts[6]));

                doZ(middle.lower, middle.upper, middle.shape());


                auto const corners
                    = for_N_make_array<8>([&](auto i) { return shift(corner, shifts[i]); });

                auto const remaining = std::apply(
                    [&](auto&&... cccs) { return box.remove_all(middle, cccs...); }, corners);

                for (auto const rem : remaining)
                    doZ(rem.lower, rem.upper, rem.shape());
            }
        }
        else
        {
            tiles.emplace_back(this->box, args...);

            // PHARE_LOG_LINE_SS(rem_per_dim);
            // for (std::uint16_t i = 0; i < tiles_per_dim[0]; ++i)
            //     PHARE_LOG_LINE_SS(starts[0][i]);
            // for (std::uint16_t i = 0; i < tiles_per_dim[1]; ++i)
            //     PHARE_LOG_LINE_SS(starts[1][i]);
            // for (std::uint16_t i = 0; i < tiles_per_dim[2]; ++i)
            //     PHARE_LOG_LINE_SS(starts[2][i]);
        }
    }
}

template<typename TileSet_t>
struct Tiler<TileSet_t, 3> : public ATiler<TileSet_t>
{
    template<typename... Args>
    auto f(Args&&... args);
};

template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t, 3>::f(Args&&... args)
{
    auto constexpr static dim          = TileSet_t::dimension;
    using Box_t                        = TileSet_t::Box_t;
    auto const& shape                  = this->box.shape();
    auto& tiles                        = this->tile_set();
    std::size_t const min_before_split = 8;

    auto split_1d = [](std::size_t const length) {
        if (length < min_before_split)
            return std::vector<std::size_t>{length};
        if (length == 8)
            return std::vector<std::size_t>{4, 4};
        if (length == 9)
            return std::vector<std::size_t>{5, 4};
        if (length == 10)
            return std::vector<std::size_t>{5, 5};
        if (length < 13)
            return std::vector<std::size_t>{6, length - 6};

        std::size_t const border = 4;
        std::size_t middle       = length - (2 * border);
        std::vector<std::size_t> sizes{border};

        if (middle / 4 < 3)
            sizes.emplace_back(middle);
        else
        {
            if (middle / 6 < 3)
            {
                sizes.emplace_back(6);
                sizes.emplace_back(middle - 6);
            }
            else
            {
                auto rem = middle % 6;
                while (middle > 5)
                {
                    std::size_t add = rem > 1 ? 2 : rem > 0 ? 1 : 0;
                    sizes.emplace_back(6 + add);
                    middle -= (6 + add);
                    rem -= add;
                }
                assert(middle == 0);
            }
        }

        sizes.push_back(border);
        return sizes;
    };

    auto const get_ranges = [](auto const& sizes) {
        std::vector<std::pair<std::size_t, std::size_t>> ranges;
        ranges.reserve(sizes.size());
        std::size_t current = 0;
        for (auto const s : sizes)
        {
            ranges.emplace_back(current, current + s);
            current += s;
        }
        return ranges;
    };

    auto const& box = this->box;

    auto const subdivide_1d = [&](auto const size_x) {
        auto const x_ranges = get_ranges(split_1d(size_x));
        tiles.reserve(x_ranges.size());
        for (auto const& [xlo, xhi] : x_ranges)
            tiles.emplace_back(Box_t{Point{xlo} + box.lower, Point{xhi} + box.lower - 1}, args...);
    };

    auto const subdivide_2d = [&](auto const size_x, auto const size_y) {
        auto const x_ranges = get_ranges(split_1d(size_x));
        auto const y_ranges = get_ranges(split_1d(size_y));
        tiles.reserve(x_ranges.size() * y_ranges.size());
        for (auto const& [xlo, xhi] : x_ranges)
            for (auto const& [ylo, yhi] : y_ranges)
                tiles.emplace_back(
                    Box_t{Point{xlo, ylo} + box.lower, Point{xhi, yhi} + box.lower - 1}, args...);
    };

    auto const subdivide_3d = [&](auto const size_x, auto const size_y, auto const size_z) {
        auto const x_ranges = get_ranges(split_1d(size_x));
        auto const y_ranges = get_ranges(split_1d(size_y));
        auto const z_ranges = get_ranges(split_1d(size_z));
        tiles.reserve(x_ranges.size() * y_ranges.size() * z_ranges.size());
        std::size_t count = 0;
        for (auto const& [xlo, xhi] : x_ranges)
            for (auto const& [ylo, yhi] : y_ranges)
                for (auto const& [zlo, zhi] : z_ranges)
                    tiles.emplace_back(Box_t{Point{xlo, ylo, zlo} + box.lower,
                                             Point{xhi, yhi, zhi} + box.lower - 1},
                                       args...),
                        ++count;

        assert(tiles.capacity() == count);
    };

    if constexpr (dim == 1)
        subdivide_1d(shape[0]);
    if constexpr (dim == 2)
        subdivide_2d(shape[0], shape[1]);
    if constexpr (dim == 3)
        subdivide_3d(shape[0], shape[1], shape[2]);

    assert(
        !any_overlaps_in(tiles, [](auto const& tile) { return static_cast<Box<int, dim>>(tile); }));
}

template<std::uint8_t impl = 3, typename TileSet_t, typename... Args>
void tile_set_make_tiles(TileSet_t& tile_set, Args&&... args)
{
    Tiler<TileSet_t, impl>{tile_set}.f(args...);
}

template<typename TileSet_t, typename TileSet0>
TileSet_t tile_set_make_from_tiles(TileSet_t const& from, auto&&... args)
{
    TileSet_t out;
    for (auto const& in : from())
        out.tiles.emplace_back(in, args...);
    return out;
}

} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP*/
