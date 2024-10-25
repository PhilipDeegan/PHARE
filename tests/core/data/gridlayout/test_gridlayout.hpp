#ifndef TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/utilities/types.hpp"


template<typename GridLayout>
class TestGridLayout : public GridLayout
{ // to expose a default constructor
public:
    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(double dl, std::uint32_t cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(dl),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    TestGridLayout(std::uint32_t cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    TestGridLayout(PHARE::core::Box<int, dim> const& amrbox)
        : GridLayout{PHARE::core::ConstArray<double, dim>(.1),
                     amrbox.shape().template toArray<std::uint32_t>(),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)},
                     amrbox}
    {
    }

    auto static make(std::uint32_t cells) { return TestGridLayout{cells}; }
    auto static make(double dl, std::uint32_t cells) { return TestGridLayout{dl, cells}; }
    auto static make(PHARE::core::Box<int, dim> const& amrbox) { return TestGridLayout{amrbox}; }

    GridLayout& operator*() { return *this; }
    GridLayout const& operator*() const { return *this; }
};

#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/
