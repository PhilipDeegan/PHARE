#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/def/phare_config.hpp"
#include "core/utilities/span.hpp"


#include <algorithm>
#include <string>
#include <tuple>
#include <cstddef>


namespace PHARE::core
{

template<typename NdArray_t>
class FieldTile : public Box<std::int32_t, NdArray_t::dimension>
{
public:
    auto constexpr static dimension = NdArray_t::dimension;
    using Super                     = Box<std::int32_t, dimension>;
    using value_type                = NdArray_t;

    template<typename... Args>
    FieldTile(Super const& box, Args&&... args)
        : Super{box}
        , data{args...}
    {
    }

    auto& operator()() { return data; }
    auto& operator()() const { return data; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

private:
    NdArray_t data;
};



template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity,
         bool hasGhosts_ = true>
class GridTile : public Box<std::int32_t, NdArray_t::dimension>
{
    // using Field_t = Field<NdArray_t::dimension, PhysicalQuantity>;

    auto static make_grid(auto&&... args)
    {
        auto const& [box, tile_layout, qty] = std::forward_as_tuple(args...);

        if constexpr (hasGhosts)
            return NdArray_t{tile_layout.allocSize(qty)};
        else
            return NdArray_t{*(box.shape().as_unsigned() - 1)};
    }


public:
    auto constexpr static hasGhosts = hasGhosts_;
    auto constexpr static dimension = NdArray_t::dimension;
    using Super                     = Box<std::int32_t, dimension>;
    using value_type                = NdArray_t;

    GridTile(Super const& box, GridLayout_t const& patch_layout, PhysicalQuantity const pq)
        : Super{box}
        , layout_{patch_layout.copy_as(box)}
        , data{make_grid(box, layout_, pq)}
    {
    }
    GridTile(GridTile&&)      = default;
    GridTile(GridTile const&) = delete;

    auto& layout() const { return layout_; }
    auto& operator()() { return data; }
    auto& operator()() const { return data; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

private:
    GridLayout_t const layout_;
    NdArray_t data;
};



template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
class GridTileSet
    : public TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>, NdArray_t::allocator_mode>
{
public:
    using Super
        = TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>, NdArray_t::allocator_mode>;
    using type                   = NdArray_t::type;
    using physical_quantity_type = PhysicalQuantity;
    // using view_t = TileSetView<FieldTile<GridLayout_t, NdArray_t, PhysicalQuantity>>

    auto constexpr static alloc_mode = NdArray_t::allocator_mode;

    template<typename Dict_t>
    GridTileSet(Dict_t const& dict, GridLayout_t const& layout, PhysicalQuantity qty)
        : Super{layout.AMRBox(), dict["tile_size"].template to<std::size_t>(), layout, qty}
        , name_{dict["name"].template to<std::string>()}
        , qty_{qty}
    {
    }

    GridTileSet(GridTileSet&&)      = default;
    GridTileSet(GridTileSet const&) = delete;

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }
    NO_DISCARD auto& name() const { return name_; }

    template<typename T>
    void fill(T const v)
    {
        for (auto& tile : *this)
            tile().fill(v);
    }

    auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell)
    {
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - tile.lower();
        // auto const lcl_tile_box = (*tile) - this->box().lower;
        // auto const cell         = lCell - lcl_tile_box.lower;
        return tile()(*cell);
    }

    auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell) const
    {
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - tile.lower();
        return tile()(*cell);
    }

    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (*this)(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (*this)(to_point<std::uint32_t>(ijk...));
    }

    auto operator*() {} // VIEW

private:
    std::string const name_;
    PhysicalQuantity const qty_;
};




template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity,
         auto alloc_mode_ = AllocatorMode::CPU>
class FieldTileSet
    : public TileSetView<FieldTile<typename NdArray_t::View>,
                         ViewSpan<FieldTile<typename NdArray_t::View>,
                                  GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>>>
{
public:
    using Super                      = TileSetView<FieldTile<typename NdArray_t::View>,
                                                   ViewSpan<FieldTile<typename NdArray_t::View>,
                                                            GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>>>;
    using type                       = NdArray_t::type;
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = alloc_mode_;

    // template<typename Dict_t>

    // TileSetView(Box_t const& box, std::array<std::size_t, dimension> const& tile_size,
    //             std::array<std::uint32_t, dimension> const& shape, Tile* tiles,
    //             std::size_t tile_nbr, Tile** cells,
    //             std::array<std::uint32_t, dimension> const& nbr_cells)

    FieldTileSet(std::string const& name, PhysicalQuantity qty)
        : Super{{}, {}, {}, nullptr, 0, nullptr, {}} // set later
        , name_{name}
        , qty_{qty}
    {
    }

    void setBuffer(std::nullptr_t ptr) _PHARE_ALL_FN_
    {
        setBuffer(static_cast<FieldTileSet*>(nullptr));
    }

    template<typename FieldLike>
    void setBuffer(FieldLike* const field) _PHARE_ALL_FN_
    {
        auto data = field ? field->data() : nullptr;
        if (data)
            super() = (*field).as([](auto&&... args) {
                auto&& [box, tile_size, shape, tiles_data, tiles_size, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                return Super{box,        tile_size,  shape,      &(*tiles_data[0]),
                             tiles_size, cells_data, cells_shape};
            });
    }

    FieldTileSet(FieldTileSet&&)      = default;
    FieldTileSet(FieldTileSet const&) = default;

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }
    NO_DISCARD auto& name() const { return name_; }

    void zero()
    {
        if (!Super::data())
            return;

        using vec_helper = PHARE::Vector<type, alloc_mode>;
        for (auto& tile : *this)
            vec_helper::fill(tile().data(), tile().size(), 0);
    }

    auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell)
    {
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - tile.lower();
        // auto const lcl_tile_box = (*tile) - this->box().lower;
        // auto const cell         = lCell - lcl_tile_box.lower;
        return tile()(*cell);
    }

    auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell) const
    {
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - tile.lower();
        return tile()(*cell);
    }

    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk)
    //     requires(sizeof...(IJK) == 0)
    // {
    //     return super()();
    // }
    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk) const
    //     requires(sizeof...(IJK) == 0)
    // {
    //     return super()();
    // }

    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (*this)(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (*this)(to_point<std::uint32_t>(ijk...));
    }


    void setData(auto* const data) _PHARE_ALL_FN_ { /*Super::setBuffer(data);*/ }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }

    auto operator*() { return *this; }       // interop atm
    auto operator*() const { return *this; } // interop atm

    void copyData(FieldTileSet fts) {}

private:
    Super const& super() { return *this; }
    Super& super() const { return *this; }



    std::string const name_;
    PhysicalQuantity const qty_;
};



} // namespace PHARE::core



#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
