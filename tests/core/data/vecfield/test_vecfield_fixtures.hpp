#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/def/phare_config.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{

// template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
// using VecField_t = VecField<Field_t<dim, alloc_mode>, HybridQuantity>;

template<typename GridLayout_, auto alloc_mde = AllocatorMode::CPU,
         auto layout_mde = LayoutMode::AoS>
using UsableVecField = UsableTensorField<GridLayout_, /*rank=*/1, alloc_mde, layout_mde>;


} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
