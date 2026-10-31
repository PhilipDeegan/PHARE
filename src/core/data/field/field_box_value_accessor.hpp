#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_VALUE_ACCESSOR_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_VALUE_ACCESSOR_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace PHARE::core
{

/** @brief Opt-in span-wise access for an argument of GridLayout::evalOn*Box / for_each_point.
 *
 * Wrapping an argument, e.g. `FieldBoxValueAccessor{Jx}`, makes for_each_point call
 * on_span(start) once per span and pass at(point) to fn instead of the argument itself:
 * the full index is computed once per span, each point is then one add along the
 * contiguous last dimension. Unwrapped arguments are passed to fn unchanged.
 *
 * Specialized for fields here, and for derivative operators (Partial, CurlComponent)
 * in core/numerics/differential_operators.hpp.
 */
template<typename T>
class FieldBoxValueAccessor;

template<typename T>
FieldBoxValueAccessor(T&) -> FieldBoxValueAccessor<T>;



/** @brief yields a reference to the field value at point, const for a const field */
template<typename Field_t>
    requires IsField<std::remove_const_t<Field_t>>
class FieldBoxValueAccessor<Field_t>
{
    using Field_nc             = std::remove_const_t<Field_t>;
    auto constexpr static dim  = Field_nc::dimension;
    auto constexpr static last = dim - 1;
    using value_type
        = std::conditional_t<std::is_const_v<Field_t>, typename Field_nc::value_type const,
                             typename Field_nc::value_type>;

    // spans run along the last index, contiguous in memory only for C ordering
    static_assert(Field_nc::is_c_ordered, "FieldBoxValueAccessor requires a C ordered field");

public:
    FieldBoxValueAccessor(Field_t& field)
        : data_{field.data()}
        , shape_{field.shape()}
    {
    }

    // unsigned wrap-around is well defined: base_ + point[last] is the flat index of point
    void on_span(auto const& start)
    {
        base_ = NdArrayViewer<dim>::idx(shape_, start) - start[last];
    }

    value_type& at(auto const& point) const { return data_[base_ + point[last]]; }

private:
    value_type* data_;
    std::array<std::uint32_t, dim> shape_;
    std::size_t base_ = 0;
};


} // namespace PHARE::core


#endif // PHARE_CORE_DATA_FIELD_FIELD_BOX_VALUE_ACCESSOR_HPP
