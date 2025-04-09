#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def/phare_config.hpp"

#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/data/tiles/tile_set.hpp"


#include <tuple>
#include <string>
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
        , arr{args...}
    {
    }

    auto data() { return arr.data(); }
    auto data() const { return arr.data(); }
    auto& operator()() { return arr; }
    auto& operator()() const { return arr; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

private:
    NdArray_t arr;
};


template<bool hasGhosts = true>
auto static grid_cells(auto&&... args)
{
    auto const& [box, tile_layout, qty] = std::forward_as_tuple(args...);

    if constexpr (hasGhosts)
        return tile_layout.allocSize(qty);
    else
        return *(box.shape().as_unsigned() - 1);
}



template<typename GridLayout_t, typename NdArray_t, typename PhysicalQuantity>
struct GridTile : public FieldTile<typename NdArray_t::View>
{
    auto constexpr static dimension = NdArray_t::dimension;
    using Super                     = FieldTile<typename NdArray_t::View>;
    using value_type                = NdArray_t;


    GridTile(auto const& box, GridLayout_t const& tile_layout, GridLayout_t const& patch_layout,
             PhysicalQuantity const& pq)
        : Super{box, nullptr, grid_cells(box, tile_layout, pq)}
        , arr{super()().shape()}
        , layout_{tile_layout}
    {
        super()().reset(arr);
    }

    GridTile(auto const& box, GridLayout_t const& patch_layout, PhysicalQuantity const& pq)
        : GridTile{box, patch_layout.copy_as(box), patch_layout, pq}
    {
    }
    GridTile(GridTile&&)      = default;
    GridTile(GridTile const&) = delete;

    auto& layout() const { return layout_; }
    auto& operator()() { return arr; }
    auto& operator()() const { return arr; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }


private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    NdArray_t arr;
    GridLayout_t layout_;
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
            super() = (*field).as([](auto&&... args) mutable {
                auto&& [box, tile_size, shape, tiles_data, tiles_size, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                return Super{box,        tile_size,       shape,      &*tiles_data[0],
                             tiles_size, &*cells_data[0], cells_shape};
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

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }


    std::string const name_;
    PhysicalQuantity const qty_;
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

    auto constexpr static alloc_mode = NdArray_t::allocator_mode;

    template<typename Dict_t>
    GridTileSet(Dict_t const& dict, GridLayout_t const& layout, PhysicalQuantity const& qty)
        : Super{layout.AMRBox(), dict["tile_size"].template to<std::size_t>(), layout, qty}
        , View{dict["name"].template to<std::string>(), qty}
        , name_{dict["name"].template to<std::string>()}
        , qty_{qty}
    {
        View::setBuffer(&super());
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

    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (**this)(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == NdArray_t::dimension)
    {
        return (**this)(to_point<std::uint32_t>(ijk...));
    }

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

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    std::string const name_;
    PhysicalQuantity const qty_;
};




} // namespace PHARE::core



#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
