#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP


#include "core/def.hpp"
#include "core/logger.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"


#include <vector>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace PHARE::core
{
template<typename D>
struct FieldBorderSumOp : public PlusEquals<D>
{
    // only work on tiles which are on the patch border
    // prevents duplicates on ghost of inner tiles
};


template<typename Field_t>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t::type>;
    using SetEqualOp = core::Equals<value_type>;
    // using layout_mode = Field_t::layout_mode;
    // auto constexpr static alloc_mode = Field_t::alloc_mode;

public:
    auto constexpr static dimension = Field_t::dimension;


    FieldBox(Field_t& field_, Box<int, dimension> const& amr_box_,
             Box<std::uint32_t, dimension> const& lcl_box_)
        : field{field_}
        , amr_box{amr_box_}
        , lcl_box{lcl_box_}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout)
        : FieldBox{field_, layout.AMRBox(), layout.ghostBoxFor(field_)}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : FieldBox{field_, layout.AMRBox(), selection}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout, Box<int, dimension> const& selection)
        : FieldBox{field_, layout.AMRBox(), layout.AMRToLocal(selection)}
    {
    }


    template<typename Operator = SetEqualOp, typename Field_t0>
    void op(FieldBox<Field_t0> const& that);



    template<typename Operator = SetEqualOp>
    void op(value_type const val);


    void append_to(std::vector<value_type>& vec) const;

    auto& offset(auto const& offset)
    {
        offset_ = offset;
        return *this;
    }

    Field_t& field;
    Box<int, dimension> amr_box;
    Box<std::uint32_t, dimension> lcl_box;
    Point<int, dimension> offset_ = ConstArray<int, dimension>();
};


template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src);


template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src)
    requires(std::is_same_v<FieldBorderSumOp<double>, Operator>)
{
    auto const pq = dst.field.physicalQuantity();
    assert(src.field.physicalQuantity() == pq);

    auto const skip = [](auto const& field, auto const& tile) {
        auto const tbox = grow(tile, 1);
        return *(tbox * field.layout().AMRBox()) == tbox;
    };
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box
        = shift(dst.field.layout().localToAMR(dst.lcl_box), src.offset_ * -1);

    for (std::size_t sidx = 0; sidx < src.field().size(); ++sidx)
    {
        auto const& src_tile = src.field()[sidx];
        if (skip(src.field, src_tile))
            continue; // skip not border tile

        auto const src_tile_ghost_box = get_box(src_tile);
        if (auto const src_overlap = src_selection_box * src_tile_ghost_box)
            for (std::size_t didx = 0; didx < dst.field().size(); ++didx)
            {
                auto& dst_tile = dst.field()[didx];
                if (skip(dst.field, dst_tile))
                    continue; // skip not border tile

                auto const dst_tile_ghost_box = shift(get_box(dst_tile), src.offset_ * -1);

                if (auto const dst_overlap = dst_selection_box * dst_tile_ghost_box)
                {
                    auto const lcl_src_box = src_tile.layout().AMRToLocal(*src_overlap);
                    auto const lcl_dst_box
                        = dst_tile.layout().AMRToLocal(shift(*dst_overlap, src.offset_));

                    assert(lcl_src_box.size() <= src_tile().size());
                    assert(lcl_dst_box.size() <= dst_tile().size());

                    auto src_it = lcl_src_box.begin();
                    auto dst_it = lcl_dst_box.begin();
                    for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end();
                         ++src_it, ++dst_it)
                        Operator{dst_tile()(*dst_it)}(src_tile()(*src_it));
                }
            }
    }
}


template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src)
{
    auto constexpr plus_equals
        = std::is_same_v<Operator, PlusEquals<typename Operator::value_type>>;

    auto const pq = dst.field.physicalQuantity();
    assert(src.field.physicalQuantity() == pq);

    auto get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };

    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box
        = shift(dst.field.layout().localToAMR(dst.lcl_box), src.offset_ * -1);

    assert(src_selection_box == dst_selection_box);

    for (auto& dst_tile : dst.field())
        if (auto const dst_overlap = dst_selection_box * shift(get_box(dst_tile), src.offset_ * -1))
            for (auto const& src_tile : src.field())
                if (auto const src_overlap = src_selection_box * get_box(src_tile))
                    if (auto const overlap = *src_overlap * *dst_overlap)
                    {
                        auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
                        auto const lcl_dst_box
                            = dst_tile.layout().AMRToLocal(shift(*overlap, src.offset_));

                        assert(lcl_src_box.size() <= src_tile().size());
                        assert(lcl_dst_box.size() <= dst_tile().size());

                        auto src_it = lcl_src_box.begin();
                        auto dst_it = lcl_dst_box.begin();
                        for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end();
                             ++src_it, ++dst_it)
                        {
                            Operator{dst_tile()(*dst_it)}(src_tile()(*src_it));
                        }
                    }
}


template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<GridTileSet<T0s...>>& dst, FieldBox<T1s...> const& src)
{
    auto const amr_selection_box = dst.field.layout().localToAMR(dst.lcl_box);
    auto const pq                = dst.field.physicalQuantity();

    auto const src_tile_ghost_box = src.amr_box;
    if (auto const src_overlap = amr_selection_box * src_tile_ghost_box)
    {
        for (std::size_t didx = 0; didx < dst.field().size(); ++didx)
        {
            auto& dst_tile                = dst.field()[didx];
            auto const dst_tile_ghost_box = dst_tile.layout().AMRGhostBoxFor(pq);
            if (auto const overlap = *src_overlap * dst_tile_ghost_box)
            {
                auto const lcl_src_box = src.lcl_box;
                auto const lcl_dst_box = dst_tile.layout().AMRToLocal(*overlap);
                auto src_it            = lcl_src_box.begin();
                auto dst_it            = lcl_dst_box.begin();
                for (; dst_it != lcl_dst_box.end(); ++src_it, ++dst_it)
                    Operator{dst_tile()(*dst_it)}(src.field(*src_it));
            }
        }
    }
}

template<typename Field_t>
template<typename Operator, typename Field_t0>
void FieldBox<Field_t>::op(FieldBox<Field_t0> const& that)
{
    operate_on_fields<Operator>(*this, that);
}


template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src)
    requires(std::is_same_v<FieldBorderSumOp<double>, Operator>)
{
    auto& dst_field = dst.field;
    auto const pq   = dst.field.physicalQuantity();
    // assert(src.field.physicalQuantity() == pq);

    auto const skip = [](auto const& field, auto const& tile) {
        auto const tbox = grow(tile, 1);
        return *(tbox * field.layout().AMRBox()) == tbox;
    };
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box
        = shift(src.field.layout().localToAMR(dst.lcl_box), src.offset_ * -1);

    for (std::size_t sidx = 0; sidx < src.field().size(); ++sidx)
    {
        auto const& src_tile = src.field()[sidx];
        if (skip(src.field, src_tile))
            continue; // skip not border tile

        auto const src_tile_ghost_box = get_box(src_tile);
        if (auto const src_overlap = src_selection_box * src_tile_ghost_box)
        {
            auto const lcl_src_box = src_tile.layout().AMRToLocal(*src_overlap);
            auto const lcl_dst_box = dst.lcl_box;
            // = dst_tile.layout().AMRToLocal(shift(*dst_overlap, src.offset_));

            assert(lcl_src_box.size() <= src_tile().size());

            auto src_it = lcl_src_box.begin();
            auto dst_it = lcl_dst_box.begin();
            for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end(); ++src_it, ++dst_it)
                Operator{dst_field(*dst_it)}(src_tile()(*src_it));
        }
    }
}

template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src)
{
    // PHARE_LOG_LINE_SS(e);
    // static_assert(false); // never called
}



template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<T0s...>& dst, FieldBox<T1s...> const& src)
{
    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
        Operator{dst.field(*dst_it)}(src.field(*src_it));
}




template<typename Operator, typename... Args>
void set_on_fields(FieldBox<FieldTileSet<Args...>>& dst, auto const val)
{
    if (dst.field().size() == 0)
        return;

    auto const pq = dst.field.physicalQuantity();
    auto const amr_selection_box
        = dst.field()[0].layout().copy_as(dst.amr_box).localToAMR(dst.lcl_box);

    for (auto& dst_tile : dst.field())
    {
        auto const dst_tile_ghost_box = dst_tile.layout().AMRGhostBoxFor(pq);
        if (auto const overlap = amr_selection_box * dst_tile_ghost_box)
        {
            auto const lcl_dst_box = dst_tile.layout().AMRToLocal(*overlap);
            for (auto dst_it = lcl_dst_box.begin(); dst_it != lcl_dst_box.end(); ++dst_it)
                Operator{dst_tile()(*dst_it)}(val);
        }
    }
}


template<typename Operator, typename... Args>
void set_on_fields(FieldBox<Args...>& dst, auto const val)
{
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++dst_it)
        Operator{dst.field(*dst_it)}(val);
}



template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::op(value_type const val)
{
    set_on_fields<Operator>(*this, val);
}

// template<typename Field_t>
// template<typename Operator>
// void FieldBox<Field_t>::op(std::vector<value_type> const& vec, std::size_t seek)
// {
//     auto dst_it = lcl_box.begin();
//     for (; dst_it != lcl_box.end(); ++seek, ++dst_it)
//         Operator{field(*dst_it)}(vec[seek]);
// }

template<typename... Args>
void append_from_fields(FieldBox<FieldTileSet<Args...> const> const& src, auto& vec)
{
    if (src.field().size() == 0)
        return;

    auto const pq = src.field.physicalQuantity();
    auto const amr_selection_box
        = src.field()[0].layout().copy_as(src.amr_box).localToAMR(src.lcl_box);

    for (auto src_tile : src.field())
    {
        auto const src_tile_ghost_box = src_tile.layout().AMRGhostBoxFor(pq);
        if (auto const overlap = amr_selection_box * src_tile_ghost_box)
        {
            auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
            for (auto src_it = lcl_src_box.begin(); src_it != lcl_src_box.end(); ++src_it)
                vec.push_back(src_tile()(*src_it));
        }
    }
}

template<typename... Args>
void append_from_fields(FieldBox<GridTileSet<Args...> const> const& src, auto& vec)
{
    append_from_fields(FieldBox<FieldTileSet<Args...> const>{src.field, src.amr_box, src.lcl_box},
                       vec);
}


template<typename... Args>
void append_from_fields(FieldBox<Args...> const& src, auto& vec)
{
    // reserve vec before use!
    auto src_it = src.lcl_box.begin();
    for (; src_it != src.lcl_box.end(); ++src_it)
        vec.push_back(src.field(*src_it));
}


template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec) const
{
    // if constexpr (is_field_tile_set_v<Field_t>)
    //     append_from_field_tiles(*this, vec);
    // else
    append_from_fields(*this, vec);
}




template<typename... T0s, typename... T1s>
void copy_fields(FieldTileSet<T0s...>& dst, FieldTileSet<T1s...> const& src)
{
    for (std::size_t idx = 0; idx < src().size(); ++idx)
        copy_fields(dst[idx](), src[idx]());
}


template<typename... T0s, auto opts>
void copy_fields(FieldTileSet<T0s...>& dst, basic::Field<opts> const& src)
{
    assert(dst().size());
    auto const layout = dst()[0].layout().copy_as(dst.box());

    for (auto& tile : dst)
        FieldBox{tile(), tile.layout(), tile.ghost_box()}.op(
            FieldBox{src, layout, tile.ghost_box()});
    // FieldBox{dst, layout}.op(FieldBox{src, layout});
}


template<typename... T0s, auto opts>
void copy_fields(basic::Field<opts>& dst, FieldTileSet<T0s...> const& src)
{
    reduce_into(src, dst);
}


template<auto opts0, auto opts1>
void copy_fields(basic::Field<opts0>& dst, basic::Field<opts1> const& src)
{
    std::memcpy(dst.data(), src.data(),
                src.size() * sizeof(typename basic::Field<opts0>::value_type));
}


template<typename Operator, typename Grid_t, typename GridTiles_t>
auto& reduce_into_(GridTiles_t const& tiles, Grid_t& grid)
{
    auto const pq      = tiles.physicalQuantity();
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };

    grid.reshape(tiles.shape());
    grid.zero();

    auto const& patch_layout = tiles[0].layout().copy_as(tiles.box());

    for (auto const& tile : tiles())
    {
        using Tile_vt           = std::decay_t<decltype(tile)>::Super;
        auto const& tile_layout = tile.layout();
        auto const& tile_box    = get_box(tile);

        FieldBox<Grid_t>{grid, patch_layout, patch_layout.AMRToLocal(tile_box)}
            .template op<Operator>(core::FieldBox<Tile_vt const>{*tile, tile_layout, tile_box});
    }

    return grid;
}


template<typename Operator, typename Grid_t, typename GridTiles_t>
auto& reduce_single_(GridTiles_t const& tiles, Grid_t& grid)
{
    auto const pq      = tiles.physicalQuantity();
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };

    grid.reshape(tiles.shape());
    grid.zero();

    auto const& patch_layout = tiles[0].layout().copy_as(tiles.box());
    auto const ghost_box     = patch_layout.AMRGhostBoxFor(pq);
    auto const field_box     = shrink(ghost_box, patch_layout.nbrGhosts());
    auto const ss0           = make_qty_nd_span_set_from(tiles, patch_layout, pq);

    for (auto const& bix : field_box)
    {
        auto const lix      = (bix - ghost_box.lower).as_unsigned();
        auto const& tile    = *ss0.cells(lix).ptr[0];
        auto const tile_gb  = tile.layout().AMRGhostBoxFor(pq);
        auto const tile_lix = (bix - tile_gb.lower).as_unsigned();

        grid(lix) = tile()(tile_lix);
    }

    for (auto const& box : ghost_box.remove(field_box))
        for (auto const& bix : box)
        {
            auto const lix      = (bix - ghost_box.lower).as_unsigned();
            auto const& tile    = *ss0.cells(lix).ptr[0];
            auto const tile_gb  = tile.layout().AMRGhostBoxFor(pq);
            auto const tile_lix = (bix - tile_gb.lower).as_unsigned();

            grid(lix) = tile()(tile_lix);
        }

    return grid;
}

template<typename Grid_t, typename GridTiles_t>
auto& reduce_single(GridTiles_t const& tiles, Grid_t& grid)
{
    if constexpr (is_field_tile_set_v<GridTiles_t>)
        reduce_single_<Equals<typename Grid_t::value_type>>(tiles, grid);
    else
        copy_fields(grid, tiles);
    return grid;
}

template<typename Tiles>
auto reduce_single(Tiles const& input)
{
    using Grid_t = Tiles::grid_type;
    Grid_t grid{input.name(), input.physicalQuantity(), input.shape()};
    reduce_single(input, grid);
    return grid;
}

template<typename Grid_t, typename GridTiles_t>
auto& reduce_into(GridTiles_t const& tiles, Grid_t& grid)
{
    if constexpr (is_field_tile_set_v<GridTiles_t>)
        reduce_into_<PlusEquals<typename Grid_t::value_type>>(tiles, grid);
    else
        copy_fields(grid, tiles);
    return grid;
}



template<typename Tiles>
auto reduce(Tiles const& input)
{
    if constexpr (!is_field_tile_set_v<Tiles>)
        return input;
    else
    {
        using Grid_t = Tiles::grid_type;
        Grid_t grid{input.name(), input.physicalQuantity(), input.shape()};
        reduce_into(input, grid);
        return grid;
    }
}



} // namespace PHARE::core


#endif
