#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/data/grid/grid.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/box/box.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include <vector>
#include <cstddef>
#include <type_traits>

namespace PHARE::core
{

template<typename Field_t>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t::type>;
    // using layout_mode = Field_t::layout_mode;
    auto constexpr static alloc_mode = Field_t::alloc_mode;

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
        : FieldBox{field_, layout.AMRBox(), layout.ghostBoxFor(field)}
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


    template<typename Operator, typename Field_t0>
    void op(FieldBox<Field_t0> const& that);

    template<typename Operator>
    void op(std::vector<value_type> const& vec, std::size_t seek = 0);

    void append_to(std::vector<value_type>& vec);


    Field_t& field;
    Box<int, dimension> amr_box;
    Box<std::uint32_t, dimension> lcl_box;
};




template<typename Operator, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<Args...>>& dst,
                       FieldBox<GridTileSet<Args...> const> const& src)
{
    auto const amr_selection_box = src.field.layout().localToAMR(src.lcl_box);

    auto const pq = dst.field.physicalQuantity();
    assert(src.field.physicalQuantity() == pq);

    for (std::size_t sidx = 0; sidx < src.field().size(); ++sidx)
    {
        auto const& src_tile          = src.field()[sidx];
        auto const src_tile_ghost_box = src_tile.layout().AMRGhostBoxFor(pq);
        if (auto const src_overlap = amr_selection_box * src_tile_ghost_box)
        {
            for (std::size_t didx = 0; didx < dst.field().size(); ++didx)
            {
                auto& dst_tile                = dst.field()[sidx];
                auto const dst_tile_ghost_box = dst_tile.layout().AMRGhostBoxFor(pq);
                if (auto const overlap = *src_overlap * dst_tile_ghost_box)
                {
                    auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
                    auto const lcl_dst_box = dst_tile.layout().AMRToLocal(*overlap);
                    auto src_it            = lcl_src_box.begin();
                    auto dst_it            = lcl_dst_box.begin();
                    for (; dst_it != lcl_dst_box.end(); ++src_it, ++dst_it)
                        Operator{dst_tile()(*dst_it)}(src_tile()(*src_it));
                }
            }
        }
    }
}




template<typename Operator, typename NdArray_>
void operate_on_fields(FieldBox<Grid<NdArray_, HybridQuantity::Scalar>>& dst,
                       FieldBox<Grid<NdArray_, HybridQuantity::Scalar> const> const& src)
{
    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
        Operator{dst.field(*dst_it)}(src.field(*src_it));
}

template<typename Field_t>
template<typename Operator, typename Field_t0>
void FieldBox<Field_t>::op(FieldBox<Field_t0> const& that)
{
    operate_on_fields<Operator>(*this, that);
}



template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::op(std::vector<value_type> const& vec, std::size_t seek)
{
    auto dst_it = lcl_box.begin();
    for (; dst_it != lcl_box.end(); ++seek, ++dst_it)
        Operator{field(*dst_it)}(vec[seek]);
}

template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec)
{
    // reserve vec before use!
    auto src_it = lcl_box.begin();
    for (; src_it != lcl_box.end(); ++src_it)
        vec.push_back(field(*src_it));
}

} // namespace PHARE::core


#endif
