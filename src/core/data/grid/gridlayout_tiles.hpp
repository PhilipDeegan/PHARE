#ifndef PHARE_CORE_DATA_GRID_GRIDLAYOUT_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRIDLAYOUT_TILES_HPP

#include "core/def.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/grid/gridlayout.hpp"

#include "core/data/grid/grid.hpp"

#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <cstddef>
#include <algorithm>


namespace PHARE::core
{

template<typename GridLayout_t>
class TileGridLayout : public GridLayout_t
{
public:
    auto constexpr static dimension = GridLayout_t::dimension;
    using Super                     = GridLayout_t;


    template<typename Field, typename Fn, typename... Args>
    void evalOnBox(Field const& field, Fn&& fn, Args&&... args) const
    {
        evalOnAnyBox(field, domainBoxFor(field), fn, args...);
    }

    template<typename Field, typename Fn, typename... Args>
    void evalOnGhostBox(Field const& field, Fn&& fn, Args&&... args) const
    {
        evalOnAnyBox(field, ghostBoxFor(field), fn, args...);
    }

    decltype(Super::AMRBox()) opBox;
};


template<typename GridLayout_t, std::size_t dim = GridLayout_t::dimension>
class GridLayoutTile : public Box<std::int32_t, GridLayout_t::dimension>
{
    template<typename... Args>
    GridLayout_t make_tile_gridlayout(Args&&... args)
    { // Finish
        auto const& [box, layout] = std::forward_as_tuple(args...);

        TileGridLayout<GridLayout_t> copy = layout.copy_as(box);
        Point<int, GridLayout_t::dimension> loShift{ConstArray<int, dim>()};
        Point<int, GridLayout_t::dimension> upShift{ConstArray<int, dim>()};

        copy.opBox.lower += loShift;
        copy.opBox.upper += upShift;

        return copy;
    }

public:
    using Super = Box<std::int32_t, GridLayout_t::dimension>;

    template<typename... Args>
    GridLayoutTile(Super const& box, Args&&... args)
        : Super{box}
        , layout{make_tile_gridlayout(box, args...)}
    {
    }

    GridLayoutTile(GridLayoutTile&&)      = default;
    GridLayoutTile(GridLayoutTile const&) = default;

    auto& operator()() { return layout; }
    auto& operator()() const { return layout; }

    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

private:
    GridLayout_t const layout;
    // Grid_t data;
};

template<typename GridLayout_t>
class GridLayoutTileSet : public TileSet<GridLayoutTile<GridLayout_t>>
{
public:
    using Super = TileSet<GridLayoutTile<GridLayout_t>>;

    template<typename Dict_t>
    GridLayoutTileSet(Dict_t const& dict, GridLayout_t const& layout)
        : Super{layout.AMRBox(), dict["tile_size"].template to<std::size_t>(),
                dict["name"].template to<std::string>(), layout}
    {
    }


private:
};


} // namespace PHARE::core


#endif // PHARE_CORE_DATA_GRID_GRIDLAYOUT_TILES_HPP
