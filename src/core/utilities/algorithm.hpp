#ifndef PHARE_ALGORITHM_HPP
#define PHARE_ALGORITHM_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/tensorfield/tensorfield.hpp"




#include <algorithm>
#include <string>

namespace PHARE
{
namespace core
{
    template<std::uint32_t lhs, std::uint32_t rhs>
    NO_DISCARD constexpr std::uint32_t max()
    {
        if constexpr (lhs < rhs)
        {
            return rhs;
        }
        else if constexpr (lhs >= rhs)
        {
            return lhs;
        }
    }


    template<typename T>
    NO_DISCARD std::string to_str(T&& t)
    {
        return t.to_str();
    }



    template<typename Container, typename ContainedT = typename Container::value_type>
    NO_DISCARD bool notIn(ContainedT& obj, Container& list)
    {
        auto sameItem = std::find_if(std::begin(list), std::end(list), [&obj](auto& currentItem) {
            return obj->name() == currentItem->name();
        });

        return sameItem == std::end(list);
    }



} // namespace core
} // namespace PHARE




namespace PHARE::core
{


template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in, FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    std::size_t const size = in().size();
    for (std::size_t i = 0; i < size; ++i)
        std::transform(in[i]().begin(), in[i]().end(), out[i]().begin(), fn);
}


template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in0, FieldTileSet<GL, ND, PQ> const& in1,
               FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    assert(in0.isUsable());
    assert(in1.isUsable());
    assert(out.isUsable());
    std::size_t const size = in0().size();
    for (std::size_t i = 0; i < size; ++i)
    {
        auto& i0 = in0[i]();
        auto& i1 = in1[i]();
        auto& o0 = out[i]();
        std::transform(in0[i]().begin(), in0[i]().end(), in1[i]().begin(), out[i]().begin(), fn);
    }
}


template<typename GL, typename ND, typename PQ>
void average(FieldTileSet<GL, ND, PQ> const& f1, FieldTileSet<GL, ND, PQ> const& f2,
             FieldTileSet<GL, ND, PQ>& avg)
{
    core::transform(f1, f2, avg, std::plus<double>());
    core::transform(avg, avg, [](double x) { return x * 0.5; });
}


template<auto opts>
void transform(basic::Field<opts> const& in, basic::Field<opts>& out, auto const fn)
{
    std::transform(in.begin(), in.end(), out.begin(), fn);
}

template<auto opts>
void transform(basic::Field<opts> const& in0, basic::Field<opts> const& in1,
               basic::Field<opts>& out, auto const fn)
{
    std::transform(in0.begin(), in0.end(), in1.begin(), out.begin(), fn);
}


template<auto opts>
void average(basic::Field<opts> const& f1, basic::Field<opts> const& f2, basic::Field<opts>& avg)
{
    core::transform(f1, f2, avg, std::plus<double>());
    core::transform(avg, avg, [](double x) { return x * 0.5; });
}

template<typename Field_t, std::size_t rank>
void average(basic::TensorField<Field_t, rank> const& vf1,
             basic::TensorField<Field_t, rank> const& vf2, basic::TensorField<Field_t, rank>& Vavg)
{
    auto constexpr static N = detail::tensor_field_dim_from_rank<rank>();

    for (std::size_t i = 0; i < N; ++i)
        average(vf1[i], vf2[i], Vavg[i]);
}


} // namespace PHARE::core

#endif
