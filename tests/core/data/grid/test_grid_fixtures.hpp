#ifndef PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode>
using Grid_t
    = Grid<NdArrayVector<dim, double, /*c_ordering=*/true, alloc_mode>, HybridQuantity::Scalar>;

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP*/
