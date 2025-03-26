#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/tiles/tile_set.hpp"


#include <string>
#include <tuple>
#include <cstddef>


namespace PHARE::core
{

template<typename GridLayout_t, typename NdArrayImpl, typename PhysicalQuantity>
class GridTile : public Box<std::int32_t, NdArrayImpl::dimension>
{
    bool constexpr static hasGhosts_ = true;
    using Grid_t                     = Grid<NdArrayImpl, PhysicalQuantity>;

    template<typename... Args>
    Grid_t static make_grid(Args&&... args)
    {
        auto const& [box, tile_layout, name, layout, qty] = std::forward_as_tuple(args...);

        if constexpr (hasGhosts)
            // return std::make_shared<Grid_t>(name, layout.copy_as(box), qty);
            return {name, tile_layout, qty};
        else
            // return std::make_shared<Grid_t>(name, qty, *(box.shape().as_unsigned() - 1));
            return {name, qty, *(box.shape().as_unsigned() - 1)};
    }

    template<typename... Args>
    auto static make_layout(Args&&... args)
    {
        auto const& [box, _, layout, __] = std::forward_as_tuple(args...);
        return layout.copy_as(box);
    }

public:
    auto constexpr static hasGhosts = hasGhosts_;
    auto constexpr static dimension = NdArrayImpl::dimension;
    using Super                     = Box<std::int32_t, dimension>;
    using value_type                = Grid_t;

    template<typename... Args>
    GridTile(Super const& box, Args&&... args)
        : Super{box}
        , layout_{make_layout(box, args...)}
        , data{make_grid(box, layout_, args...)}
    {
    }

    GridTile(GridTile&&)      = default;
    GridTile(GridTile const&) = delete;

    auto& operator()() { return data; }
    auto& operator()() const { return data; }

    auto& layout() const { return layout_; }

    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

private:
    GridLayout_t const layout_;
    Grid_t data;
};

template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
class GridTileSet
    : public TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>, NdArray_t::allocator_mode>
{
public:
    using Super
        = TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>, NdArray_t::allocator_mode>;

    auto constexpr static alloc_mode = NdArray_t::allocator_mode;

    template<typename Dict_t>
    GridTileSet(Dict_t const& dict, GridLayout_t const& layout, PhysicalQuantity qty)
        : Super{layout.AMRBox(), dict["tile_size"].template to<std::size_t>(),
                dict["name"].template to<std::string>(), layout, qty}
        , name_{dict["name"].template to<std::string>()}
        , qty_{qty}
    {
    }

    GridTileSet(GridTileSet&&)      = default;
    GridTileSet(GridTileSet const&) = delete;

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }

    template<typename T>
    void fill(T const v)
    {
        for (auto& tile : *this)
            tile().fill(v);
    }

    // auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell)
    // {
    //     PHARE_LOG_LINE_SS(name_ << " " << lCell);

    //     assert(Super::at(lCell));
    //     auto& tile = *Super::at(lCell);

    //     // PHARE_LOG_LINE_SS(this->box());
    //     // PHARE_LOG_LINE_SS((*tile));

    //     auto const lcl_tile_box = (*tile) - this->box().lower;
    //     // PHARE_LOG_LINE_SS(lcl_tile_box);
    //     // PHARE_LOG_LINE_SS(tile.lower);

    //     auto const cell = lCell - lcl_tile_box.lower;
    //     PHARE_LOG_LINE_SS(cell);
    //     return tile()(*cell);
    // }

    // auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell) const
    // {
    //     assert(Super::at(lCell));
    //     auto& tile      = *Super::at(lCell);
    //     auto const cell = lCell - tile.lower();
    //     return tile()(*cell);
    // }

    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk)
    // {
    //     return (*this)(to_point<std::uint32_t>(ijk...));
    // }
    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk) const
    // {
    //     return (*this)(to_point<std::uint32_t>(ijk...));
    // }

private:
    std::string const name_;
    PhysicalQuantity const qty_;
};


} // namespace PHARE::core


#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
