#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP

#include "core/data/field/field.hpp"

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using Field_t = Field<dim, HybridQuantity::Scalar, double, alloc_mode>;

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP*/
