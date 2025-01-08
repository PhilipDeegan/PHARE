#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "core/data/vecfield/vecfield.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using VecField_t = VecField<Field_t<dim, alloc_mode>, HybridQuantity>;

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using UsableVecField = UsableTensorField<dim, /*rank=*/1, alloc_mode>;


} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
