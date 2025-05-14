#ifndef PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"

namespace PHARE::core
{


template<bool binary_eq = false, typename GridLayout_, typename NdArray0, typename NdArray1>
EqualityReport compare_fields(GridTileSet<GridLayout_, NdArray0, HybridQuantity::Scalar> const& ref,
                              Grid<NdArray1, HybridQuantity::Scalar> const& cmp,
                              double const diff = 1e-15)
{
    auto const same_sizes = ref.size() == cmp.size();

    if (!same_sizes)
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};

    auto tmp = cmp;
    tmp.zero();
    reduce_into(ref, tmp);
    return compare_fields(*tmp, *cmp, diff);
}


template<bool binary_eq = false, typename GridLayout_, typename NdArray0, typename NdArray1>
EqualityReport compare_fields(Grid<NdArray0, HybridQuantity::Scalar> const& ref,
                              GridTileSet<GridLayout_, NdArray1, HybridQuantity::Scalar> const& cmp,
                              double const diff = 1e-15)
{
    auto const same_sizes = ref.size() == cmp.size();

    if (!same_sizes)
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};

    auto tmp = ref;
    tmp.zero();
    reduce_into(cmp, tmp);
    return compare_fields(*ref, *tmp, diff);
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP*/
