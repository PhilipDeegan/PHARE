
#include "phare_core.hpp"
#include "phare/phare.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "amr/data/particles/particles_data.hpp"


#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "tests/amr/test_hierarchy_fixtures.hpp"

#include "gtest/gtest.h"

#include <SAMRAI/pdat/CellGeometry.h>

static constexpr std::size_t interp = 1;

template<std::size_t _dim>
struct TestParam
{
    auto constexpr static dim = _dim;
    using PhareTypes          = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Box_t               = PHARE::core::Box<int, dim>;
    using ParticleArray_t     = PHARE::core::ParticleArray<dim>;

    using Hierarchy_t
        = AfullHybridBasicHierarchy<dim, interp, defaultNbrRefinedParts<dim, interp>()>;
};

template<typename TestParam>
struct ParticleScheduleHierarchyTest : public ::testing::Test
{
};

using ParticlesDatas = testing::Types<TestParam<1>, TestParam<2> /*,TestParam<3>*/>;


TYPED_TEST_SUITE(ParticleScheduleHierarchyTest, ParticlesDatas);
namespace PHARE::core
{


TYPED_TEST(ParticleScheduleHierarchyTest, testing)
{
    //

    //
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
