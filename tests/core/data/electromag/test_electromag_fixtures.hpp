#ifndef PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP

#include "phare_core.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include <cmath>
#include <cassert>
#include <functional>

namespace PHARE::core
{

template<typename EM, typename GridLayout>
auto default_em_init(EM& em, GridLayout const& layout)
{
    auto constexpr static dim = GridLayout::dimension;

    auto setter = [&](auto& v) {
        auto box = layout.ghostBoxFor(v);
        for (std::size_t i = 0; i < box.upper[0]; i += 2)
        {
            if constexpr (dim == 1)
            {
                v(i) += .05;
            }
            else
            {
                for (std::size_t j = 0; j < box.upper[1]; j += 2)
                {
                    if constexpr (dim == 2)
                    {
                        v(i, j) += .05;
                    }
                    else
                    {
                        if constexpr (dim == 3)
                        {
                            for (std::size_t k = 0; k < box.upper[2]; k += 2)
                            {
                                v(i, j, k) += .05;
                            }
                        }
                    }
                }
            }
        }
    };
    for (auto& xyz : em.E)
        setter(*xyz);
    for (auto& xyz : em.B)
        setter(*xyz);
}


template<std::size_t dim>
class UsableElectromag : public Electromag<VecField_t<dim>>
{
    void _set()
    {
        E.set_on(Super::E);
        B.set_on(Super::B);
    }

public:
    using Super = Electromag<VecField_t<dim>>;

    template<typename GridLayout>
    UsableElectromag(GridLayout const& layout)
        : Super{"EM"}
        , E{"EM_E", layout, HybridQuantity::Vector::E, .1}
        , B{"EM_B", layout, HybridQuantity::Vector::B, .1}
    {
        default_em_init(*this, layout);
        _set();
    }

    UsableElectromag(UsableElectromag&& that)
        : Super{std::forward<Super>(that)}
        , E{std::move(that.E)}
        , B{std::move(that.B)}
    {
        _set();
    }


    Super& view() { return *this; }
    Super const& view() const { return *this; }
    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }


    UsableVecField<dim> E, B;
};


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP */
