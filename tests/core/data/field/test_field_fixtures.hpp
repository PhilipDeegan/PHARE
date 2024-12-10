#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP

#include "core/data/field/field.hpp"
#include "core/utilities/equality.hpp"

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using Field_t = Field<dim, HybridQuantity::Scalar, double, alloc_mode>;


template<bool binary_eq = false, std::size_t dim, auto am0, auto am1>
EqualityReport compare_fields(Field<dim, HybridQuantity::Scalar, double, am0> const& ref,
                              Field<dim, HybridQuantity::Scalar, double, am1> const& cmp,
                              [[maybe_unused]] double const diff)
{
    auto const float_eq = [&](auto const a, auto const b) {
        if constexpr (binary_eq)
            return a == b;
        else
            return float_equals(a, b, diff);
    };

    std::stringstream log;

    auto const& ref_dat = ref.data();
    auto const& cmp_dat = cmp.data();
    std::size_t eqvals = 0, eqnot0 = 0;
    for (std::size_t i = 0; i < ref.size(); ++i)
        if (float_eq(ref_dat[i], cmp_dat[i]))
        {
            ++eqvals;
            if (ref_dat[i] != 0 and cmp_dat[i] != 0)
                ++eqnot0;
        }
    if (eqvals != ref.size())
    {
        auto const bad = ref.size() - eqvals;
        log << "Field value mismatch: \n";
        log << "ok(" << eqvals << ") - ";
        log << "ok!=0(" << eqnot0 << ") - ";
        log << "bad(" << bad << ")\n";
        return EqualityReport{false, log.str()};
    }



    return EqualityReport{true};
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP*/
