#ifndef PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/field/field_box.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"

namespace PHARE::core
{


template<typename... T0s, typename... T1s>
EqualityReport compare_fields(GridTileSet<T0s...> const& ref, Grid<T1s...> const& cmp,
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


template<typename... T0s, typename... T1s>
EqualityReport compare_fields(Grid<T0s...> const& ref, GridTileSet<T1s...> const& cmp,
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
