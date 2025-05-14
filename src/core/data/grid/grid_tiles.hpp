#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def/phare_config.hpp"


#include "core/def.hpp"
#include "core/data/grid/grid.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"


#include <tuple>
#include <string>
#include <cstddef>
#include <stdexcept>


namespace PHARE::core
{

template<typename Field_t>
class FieldBox;

template<typename GridLayout_t, typename Field_t>
class FieldTile : public Box<std::int32_t, GridLayout_t::dimension>
{
public:
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = Field_t::alloc_mode;
    using Super                      = Box<std::int32_t, dimension>;
    using value_type                 = Field_t;
    using type                       = Field_t::type;

    template<typename... Args>
    FieldTile(GridLayout_t const& layout, Field_t const& field)
        : Super{layout.AMRBox()}
        , field_{field}
        , layout_{layout}
    {
        // PHARE_LOG_LINE_SS(layout.AMRBox());
    }

    auto& operator()() _PHARE_ALL_FN_ { return field_; }
    auto& operator()() const _PHARE_ALL_FN_ { return field_; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }
    auto& layout() const _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_SS(layout_.AMRBox());
        return layout_;
    }


    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) _PHARE_ALL_FN_
    {
        return field_(point);
    }
    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) const _PHARE_ALL_FN_
    {
        return field_(point);
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == dimension)
    _PHARE_ALL_FN_
    {
        return field_(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == dimension)
    _PHARE_ALL_FN_
    {
        return field_(to_point<std::uint32_t>(ijk...));
    }

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return field_.physicalQuantity(); }

    bool isUsable() const { return field_.isUsable(); }
    bool isSettable() const { return !isUsable(); }

private:
    Field_t field_;
    GridLayout_t layout_;
};



template<bool hasGhosts = true>
auto static grid_cells(auto&&... args)
{
    auto const& [layout, qty] = std::forward_as_tuple(args...);

    if constexpr (hasGhosts)
        return layout.allocSize(qty);
    else
        return *(layout.AMRBox().shape().as_unsigned() - 1);
}



template<typename GridLayout_t, typename Grid_t, typename Field_t>
struct GridTile : public FieldTile<GridLayout_t, Field_t>
{
    auto constexpr static dimension = GridLayout_t::dimension;
    using Super                     = FieldTile<GridLayout_t, Field_t>;
    using value_type                = Grid_t::Super;
    using NdArray_t                 = Grid_t::Super;
    using Super::operator();

    GridTile(auto const& layout, auto const& pq)
        : Super{layout, Field_t{pq}}
        , arr{grid_cells(layout, pq)}
    {
        reset();
        // PHARE_LOG_LINE_SS(layout.AMRBox());
        assert(!Super::layout().AMRBox().isEmpty());
    }

    GridTile(auto const& box, auto const& patch_layout, auto const& pq) // called when making tiles
        : GridTile(patch_layout.copy_as(box), pq)
    {
        // PHARE_LOG_LINE_SS(box);
    }

    GridTile(GridTile const& that)
        : Super{that.layout(), Field_t{(*that).physicalQuantity()}}
        , arr{that.arr}
    {
        reset();
    }

    GridTile(GridTile&&)                 = delete;
    GridTile& operator=(GridTile const&) = delete;
    GridTile& operator=(GridTile&&)      = delete;


    // auto& operator()() { return arr; }
    // auto& operator()() const { return arr; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return (**this).physicalQuantity(); }


private:
    void reset() { (**this)() = Field_t{(**this).physicalQuantity(), arr.data(), arr.shape()}; }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

    NdArray_t arr;
};

template<typename GridLayout_t, typename Grid_t, typename Field_t>
struct FieldTileSetter
{
    using Tile_t     = GridTile<GridLayout_t, Grid_t, Field_t>;
    using Tile_vt    = FieldTile<GridLayout_t, Field_t>;
    using Span_t     = ViewSpan<Tile_vt, Tile_t>;
    using Span_pt    = ViewSpan<Tile_vt*, Tile_t*>;
    using value_type = TileSetView<Tile_vt, Span_t, Span_pt>;
};

template<typename GridLayout_t, typename Grid_t, typename Field_t>
class FieldTileSet : public FieldTileSetter<GridLayout_t, Grid_t, Field_t>::value_type
{
    using local_types = FieldTileSetter<GridLayout_t, Grid_t, Field_t>;

public:
    using grid_type  = Grid_t;
    using field_type = Field_t;
    using Super      = local_types::value_type;
    using value_type = local_types::Tile_vt;

    using physical_quantity_type     = Grid_t::physical_quantity_type;
    using type                       = Grid_t::type;
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = Grid_t::alloc_mode;

    FieldTileSet(std::string const& name, auto physicalQuantity)
        : Super{{}, nullptr, 0, nullptr, {}} // set later
        , name_{name}
        , qty_{physicalQuantity}
    {
    }

    FieldTileSet(FieldTileSet const&)            = default;
    FieldTileSet(FieldTileSet&&)                 = delete;
    FieldTileSet& operator=(FieldTileSet const&) = delete;
    FieldTileSet& operator=(FieldTileSet&&)      = delete;


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
            // assert((**this)[0]().data());
            super() = field->super().as([](auto&&... args) mutable {
                auto&& [box, tiles_data, n_tiles, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                // PHARE_LOG_LINE_SS(n_tiles);
                assert(cells_data);
                return Super{box, &*tiles_data[0], n_tiles,
                             reinterpret_cast<FieldTile<GridLayout_t, Field_t>**>(cells_data),
                             cells_shape};
            });
            check();
            ghost_box_ = field->layout().AMRGhostBoxFor(qty_);
        }
        else
            super() = Super{{}, nullptr, 0, nullptr, {}};
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

    // auto& operator()(Point<std::uint32_t, dimension> const& lCell)
    // {
    //     assert(this->isUsable());
    //     assert(Super::at(lCell));
    //     auto& tile      = *Super::at(lCell);
    //     auto const cell = lCell - (tile.lower - Super::box().lower);
    //     return tile()(*cell);
    // }

    // auto& operator()(Point<std::uint32_t, dimension> const& lCell) const
    // {
    //     assert(this->isUsable());
    //     assert(Super::at(lCell));
    //     auto& tile      = *Super::at(lCell);
    //     auto const cell = lCell - (tile.lower - Super::box().lower);
    //     return tile()(*cell);
    // }

    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk)
    //     requires(sizeof...(IJK) == dimension)
    // {
    //     return (*this)(to_point<std::uint32_t>(ijk...));
    // }
    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk) const
    //     requires(sizeof...(IJK) == dimension)
    // {
    //     return (*this)(to_point<std::uint32_t>(ijk...));
    // }

    auto& operator()()
    {
        assert(isUsable());
        return super()();
    }
    auto& operator()() const
    {
        assert(isUsable());
        return super()();
    }


    void setData(auto* const data) _PHARE_ALL_FN_
    {
        throw std::runtime_error("fix?");
        /*Super::setBuffer(data);*/
    }

    bool isUsable() const
    {
        auto const b = Super::data() != nullptr;
        if (b)
            check();
        return b;
    }
    bool isSettable() const { return !isUsable(); }

    auto& operator*() { return *this; }       // interop atm
    auto& operator*() const { return *this; } // interop atm

    void copyData(FieldTileSet const& that)
    {
        for (std::size_t tidx = 0; tidx < super().size(); ++tidx)
            std::copy(that[tidx]().data(), that[tidx]().data() + that[tidx]().size(),
                      super()[tidx]().data());
    }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return ghost_box_.size(); } // NOT ntiles!
    NO_DISCARD auto shape() const _PHARE_ALL_FN_ { return *ghost_box_.shape().as_unsigned(); }
    NO_DISCARD auto ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }


    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }
    // NO_DISCARD auto size() const _PHARE_ALL_FN_ { return layout_.allocSize(qty_); }

private:
    void check() const
    {
        for (auto& tile : super())
        {
            assert(tile().size() < static_cast<std::size_t>(1e6));
        }
    }

    Super& super() { return *this; }
    Super const& super() const { return *this; }


    std::string name_;
    physical_quantity_type qty_;
    Box<int, dimension> ghost_box_{};
};


template<typename GridLayout_t, typename Grid_t, typename Field_t>
class GridTileSet : public TileSet<GridTile<GridLayout_t, Grid_t, Field_t>, Grid_t::alloc_mode>,
                    public FieldTileSet<GridLayout_t, Grid_t, Field_t>

{
public:
    using grid_type  = Grid_t;
    using field_type = Field_t;
    using Super      = TileSet<GridTile<GridLayout_t, Grid_t, Field_t>, Grid_t::alloc_mode>;
    using View       = FieldTileSet<GridLayout_t, Grid_t, Field_t>;
    using type       = Grid_t::value_type;
    using physical_quantity_type = Grid_t::physical_quantity_type;
    using value_type             = GridTile<GridLayout_t, Grid_t, Field_t>;
    using Super::operator[];
    using View::operator();
    using View::box;

    auto constexpr static alloc_mode = Grid_t::alloc_mode;
    auto constexpr static dimension  = Grid_t::dimension;

    GridTileSet(std::string const& name, GridLayout_t const& layout,
                physical_quantity_type const qty)
        : Super{layout.AMRBox(), layout, qty}
        , View{name, qty}
        , layout_{layout}
    {
        assert(this->super()[0]().data());
        View::setBuffer(this);
        assert(View::isUsable());
        assert((**this)[0]().data());
        // PHARE_LOG_LINE_SS(layout_.AMRBox());
    }

    GridTileSet(GridTileSet const& that)
        : Super{that.super()}
        , View{that}
        , layout_{that.layout_}
    {
        View::setBuffer(this);
        assert((**this)[0]().data());
        assert(View::isUsable());
    }

    GridTileSet(GridTileSet&&) = delete;


    // FieldTileSet(FieldTileSet const&)            = delete;
    // FieldTileSet(FieldTileSet&&)                 = delete;
    GridTileSet& operator=(GridTileSet const&) = delete;
    GridTileSet& operator=(GridTileSet&&)      = delete;

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
    //     requires(sizeof...(IJK) == dimension)
    // {
    //     return (**this)(to_point<std::uint32_t>(ijk...));
    // }
    // template<typename... IJK>
    // auto& operator()(IJK const&... ijk) const
    //     requires(sizeof...(IJK) == dimension)
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
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return View::size(); }
    auto& layout() const { return layout_; }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

private:
    GridLayout_t layout_;
};


template<typename Grid_t, typename GridTiles_t>
auto& reduce_into(GridTiles_t const& tiles, Grid_t& grid)
{
    // if (tiles.size() > grid.vector().capacity())
    //     throw std::runtime_error("Cannot reduce into grid of smaller size");

    grid.reshape(tiles.shape());
    grid.zero();

    auto const& patch_layout = tiles[0].layout().copy_as(tiles.box());

    for (auto const& tile : tiles())
    {
        using Tile_vt           = std::decay_t<decltype(tile)>::Super;
        auto const& tile_layout = tile.layout();
        auto const& tile_box    = shrink(tile_layout.AMRGhostBoxFor(tiles.physicalQuantity()),
                                         patch_layout.nbrGhosts());

        FieldBox<Grid_t>{grid, patch_layout, tile_box}.template op</*Plus*/ Equals<double>>(
            core::FieldBox<Tile_vt const>{*tile, tile_layout, tile_box});
    }

    return grid;
}


template<typename T>
struct is_field_tile_set : std::false_type
{
};


template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<FieldTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<GridTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename T>
auto static constexpr is_field_tile_set_v = is_field_tile_set<T>::value;




template<typename Tiles>
auto reduce(Tiles const& input)
{
    if constexpr (!is_field_tile_set_v<Tiles>)
        return input;
    else
    {
        using Grid_t = Tiles::grid_type;
        // using PhysicalQuantity = Tiles::physical_quantity_type;
        Grid_t grid{input.name(), input.physicalQuantity(), input.shape()};
        reduce_into(input, grid);
        return grid;
    }
}


} // namespace PHARE::core



#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
