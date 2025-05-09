#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def/phare_config.hpp"

#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
// #include "core/data/field/field_box.hpp"
#include "core/utilities/types.hpp"


#include <tuple>
#include <string>
#include <cstddef>


namespace PHARE::core
{

template<typename Field_t>
class FieldBox;

template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity, auto alloc_mde>
class FieldTile : public Box<std::int32_t, NdArray_t::dimension>
{
public:
    auto constexpr static dimension  = NdArray_t::dimension;
    auto constexpr static alloc_mode = alloc_mde;
    using Super                      = Box<std::int32_t, dimension>;
    using value_type                 = NdArray_t;
    using type                       = NdArray_t::type;

    template<typename... Args>
    FieldTile(Super const& box, GridLayout_t const& layout, auto ptr, auto shape,
              PhysicalQuantity const& pq)
        : Super{box}
        , arr{ptr, shape}
        , layout_{layout}
        , qty_{pq}
    {
        // PHARE_LOG_LINE_SS(layout_.AMRBox());
    }

    auto data() _PHARE_ALL_FN_ { return arr.data(); }
    auto data() const _PHARE_ALL_FN_ { return arr.data(); }
    auto size() const _PHARE_ALL_FN_ { return arr.size(); } // *deref* to get box size
    auto& operator()() _PHARE_ALL_FN_ { return arr; }
    auto& operator()() const _PHARE_ALL_FN_ { return arr; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }
    auto& layout() const _PHARE_ALL_FN_ { return layout_; }
    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }



    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) _PHARE_ALL_FN_
    {
        return arr(point);
    }
    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) const _PHARE_ALL_FN_
    {
        return arr(point);
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == NdArray_t::dimension)
    _PHARE_ALL_FN_
    {
        return arr(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == NdArray_t::dimension)
    _PHARE_ALL_FN_
    {
        return arr(to_point<std::uint32_t>(ijk...));
    }

private:
    NdArray_t arr;
    GridLayout_t layout_;
    PhysicalQuantity qty_;
};



template<bool hasGhosts = true>
auto static grid_cells(auto&&... args)
{
    auto const& [box, tile_layout, qty] = std::forward_as_tuple(args...);

    if constexpr (hasGhosts)
    {
        return tile_layout.allocSize(qty);
    }
    else
        return *(box.shape().as_unsigned() - 1);
}



template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
struct GridTile : public FieldTile<GridLayout_t, typename NdArray_t::View, PhysicalQuantity,
                                   NdArray_t::allocator_mode>
{
    auto constexpr static dimension = NdArray_t::dimension;
    using Super      = FieldTile<GridLayout_t, typename NdArray_t::View, PhysicalQuantity,
                                 NdArray_t::allocator_mode>;
    using value_type = NdArray_t;


    GridTile(auto const& box, GridLayout_t const& tile_layout, GridLayout_t const& patch_layout,
             PhysicalQuantity const& pq)
        : Super{box, tile_layout, nullptr, grid_cells(box, tile_layout, pq), pq}
        , arr{super()().shape()}
    {
        super()().reset(arr);
        assert(!Super::layout().AMRBox().isEmpty());
    }

    GridTile(auto const& box, GridLayout_t const& patch_layout, PhysicalQuantity const& pq)
        : GridTile{box, patch_layout.copy_as(box), patch_layout, pq}
    {
        // PHARE_LOG_LINE_SS(box);
        // PHARE_LOG_LINE_SS(patch_layout.AMRBox());
    }
    GridTile(GridTile const& that)
        : Super{that}
        , arr{super()().shape()}
    {
        super()().reset(arr);
        assert(!Super::layout().AMRBox().isEmpty());
    }

    GridTile(GridTile&&) = default;


    auto& operator()() { return arr; }
    auto& operator()() const { return arr; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }



    auto& operator()(std::array<std::uint32_t, dimension> const& idx) { return arr(idx); }
    auto& operator()(std::array<std::uint32_t, dimension> const& idx) const { return arr(idx); }


    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return arr(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return arr(to_point<std::uint32_t>(ijk...));
    }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    NdArray_t arr;
};

template<typename GridLayout_t, typename NdArray_t, typename NdArray_vt, typename PhysicalQuantity>
struct FieldTileSetter
{
    using Tile_t = GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>;
    using Tile_vt
        = FieldTile<GridLayout_t, NdArray_vt, PhysicalQuantity, NdArray_t::allocator_mode>;
    using Span_t     = ViewSpan<Tile_vt, Tile_t>;
    using Span_pt    = ViewSpan<Tile_vt*, Tile_t*>;
    using value_type = TileSetView<Tile_vt, Span_t, Span_pt>;
};

template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
class FieldTileSet : public FieldTileSetter<GridLayout_t, NdArray_t, typename NdArray_t::View,
                                            PhysicalQuantity>::value_type
{
    using local_types
        = FieldTileSetter<GridLayout_t, NdArray_t, typename NdArray_t::View, PhysicalQuantity>;
    using NdArray_vt = NdArray_t::View;

public:
    using Super                      = local_types::value_type;
    using value_type                 = local_types::Tile_vt;
    using type                       = NdArray_t::type;
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = NdArray_t::allocator_mode;

    FieldTileSet(std::string const& name, PhysicalQuantity qty)
        : Super{{}, nullptr, 0, nullptr, {}} // set later
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
        {
            super() = (*field).as([](auto&&... args) mutable {
                auto&& [box, tiles_data, tiles_size, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                assert(cells_data);
                return Super{box, &*tiles_data[0], tiles_size,
                             reinterpret_cast<FieldTile<GridLayout_t, NdArray_vt, PhysicalQuantity,
                                                        alloc_mode>**>(cells_data),
                             cells_shape};
            });
        }
    }


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
        assert(this->isUsable());
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - (tile.lower - Super::box().lower);
        return tile()(*cell);
    }

    auto& operator()(Point<std::uint32_t, NdArray_t::dimension> const& lCell) const
    {
        assert(this->isUsable());
        assert(Super::at(lCell));
        auto& tile      = *Super::at(lCell);
        auto const cell = lCell - (tile.lower - Super::box().lower);
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

    auto& operator()() { return super()(); }
    auto& operator()() const { return super()(); }


    void setData(auto* const data) _PHARE_ALL_FN_ { /*Super::setBuffer(data);*/ }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }

    auto& operator*() { return *this; }       // interop atm
    auto& operator*() const { return *this; } // interop atm

    void copyData(FieldTileSet fts) {}

    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }
    // NO_DISCARD auto size() const _PHARE_ALL_FN_ { return layout_.allocSize(qty_); }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }


    std::string name_;
    PhysicalQuantity qty_;
};


template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
class GridTileSet : public TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>,
                                   NdArray_t::allocator_mode>,
                    public FieldTileSet<GridLayout_t, NdArray_t, PhysicalQuantity>

{
public:
    using Super
        = TileSet<GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>, NdArray_t::allocator_mode>;
    using View                   = FieldTileSet<GridLayout_t, NdArray_t, PhysicalQuantity>;
    using type                   = NdArray_t::type;
    using physical_quantity_type = PhysicalQuantity;
    using value_type             = GridTile<GridLayout_t, NdArray_t, PhysicalQuantity>;
    using Super::operator[];
    using View::operator();

    auto constexpr static alloc_mode = NdArray_t::allocator_mode;
    auto constexpr static dimension  = View::dimension;

    GridTileSet(std::string const& name, GridLayout_t const& layout, PhysicalQuantity const& qty)
        : Super{layout.AMRBox(), layout, qty}
        , View{name, qty}
        , layout_{layout}
    {
        View::setBuffer(&super());
        // PHARE_LOG_LINE_SS(layout_.AMRBox());
    }

    GridTileSet(GridTileSet const& that)
        : Super{that.super()}
        , View{that}
        , layout_{that.layout_}
    {
        View::setBuffer(&super());
    }

    GridTileSet(GridTileSet&&) = default;

    // NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return View::qty_; }
    // NO_DISCARD auto& name() const { return View::name_; }

    template<typename T>
    void fill(T const v)
    {
        for (auto& tile : *this)
            tile().fill(v);
    }

    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk)
    //     requires(sizeof...(IJK) == NdArray_t::dimension)
    // {
    //     return (**this)(to_point<std::uint32_t>(ijk...));
    // }
    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk) const
    //     requires(sizeof...(IJK) == NdArray_t::dimension)
    // {
    //     return (**this)(to_point<std::uint32_t>(ijk...));
    // }

    auto& operator()() { return super()(); }
    auto& operator()() const { return super()(); }

    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }

    View& operator*() { return *this; }
    View const& operator*() const { return *this; }

    auto data() { return super().data(); }
    auto data() const { return super().data(); }
    NO_DISCARD auto begin() { return super().begin(); }
    NO_DISCARD auto begin() const { return super().begin(); }
    NO_DISCARD auto end() { return super().end(); }
    NO_DISCARD auto end() const { return super().end(); }
    NO_DISCARD auto size() const _PHARE_ALL_FN_
    {
        return product(layout_.allocSize(View::physicalQuantity()));
    }
    auto& layout() const { return layout_; }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    GridLayout_t layout_;
};


template<typename Grid_t, typename GridTiles_t>
void reduce_into(GridTiles_t const& tiles, Grid_t& grid)
{
    if (tiles.size() != grid.size())
        throw std::runtime_error("Cannot reduce into grid of different size");

    auto const& patch_layout = tiles.layout();

    for (auto const& tile : tiles())
    {
        using Tile_vt           = std::decay_t<decltype(tile)>::Super;
        auto const& tile_layout = tile.layout();

        FieldBox<Grid_t>{grid, patch_layout, tile_layout.AMRGhostBoxFor(tiles.physicalQuantity())}
            .template op<PlusEquals<double>>(core::FieldBox<Tile_vt const>{*tile, tile_layout});
    }
}


template<typename T>
struct is_field_tile_set : std::false_type
{
};


template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<FieldTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename T>
auto static constexpr is_field_tile_set_v = is_field_tile_set<T>::value;

} // namespace PHARE::core



#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
