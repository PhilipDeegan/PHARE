#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP

#include "core/data/field/field.hpp"
#include "core/utilities/equality.hpp"

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using Field_t = Field<dim, HybridQuantity::Scalar, double, alloc_mode>;


template<bool binary_eq = false>
struct FieldComparator
{
    auto float_eq(auto const a, auto const b) const
    {
        if constexpr (binary_eq)
            return a == b;
        else
            return float_equals(a, b, diff);
    };

    template<typename F0, typename F1>
    auto operator()(F0 const& ref, F1 const& cmp)
    {
        auto const& ref_dat = ref.data();
        auto const& cmp_dat = cmp.data();
        for (std::size_t i = 0; i < ref.size(); ++i)
        {
            if (ref_dat[i] == 0)
                ++ref0;
            if (cmp_dat[i] == 0)
                ++cmp0;

            if (float_eq(ref_dat[i], cmp_dat[i]))
            {
                ++eqvals;
                if (ref_dat[i] != 0 and cmp_dat[i] != 0)
                    ++eqnot0;
            }
        }
        return std::make_tuple(eqvals, eqnot0, ref0, cmp0);
    }

    double const diff  = 1e-15;
    std::size_t eqvals = 0, eqnot0 = 0, ref0 = 0, cmp0 = 0;
};

using FloatFieldComparator_t = FieldComparator<false>;

template<bool binary_eq = false, std::size_t dim, auto am0, auto am1>
EqualityReport compare_fields(Field<dim, HybridQuantity::Scalar, double, am0> const& ref,
                              Field<dim, HybridQuantity::Scalar, double, am1> const& cmp,
                              [[maybe_unused]] double const diff)
{
    auto const same_sizes = ref.size() == cmp.size();

    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    std::stringstream log;

    auto const [eqvals, eqnot0, ref0, cmp0] = FloatFieldComparator_t{diff}(ref, cmp);

    if (eqvals != ref.size())
    {
        auto const bad = ref.size() - eqvals;
        log << "Field value mismatch: \n";
        log << "ok(" << eqvals << ") - ";
        log << "ok!=0(" << eqnot0 << ") - ";
        log << "bad(" << bad << ") - ";
        log << "ref0(" << ref0 << ") - ";
        log << "cmp0(" << cmp0 << ")\n";
        return EqualityReport{false, log.str()};
    }

    return EqualityReport{true};
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP*/
