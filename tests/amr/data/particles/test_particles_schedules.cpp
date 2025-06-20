
#include "phare_core.hpp"


#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>
#include <amr/utilities/box/amr_box.hpp>

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"



#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>

#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 100;

template<auto opts>
struct TestParam
{
    auto constexpr static dim = opts.dimension;
    using PhareTypes          = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Hierarchy_t         = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct ParticleScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    auto constexpr static dim = TestParam::dim;
    using Hierarchy_t         = TestParam::Hierarchy_t;
    using GridLayout_t        = TestParam::GridLayout_t;
    using ResourceManager_t   = Hierarchy_t::ResourcesManagerT;

    std::string configFile = "test_particles_schedules_inputs/" + std::to_string(dim) + "d_L0.txt";
    Hierarchy_t hierarchy{configFile};
};

// clang-format off
using ParticlesDatas = testing::Types<
   TestParam<SimOpts{}>

PHARE_WITH_MKN_GPU(
  ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
)

>;
// clang-format on


TYPED_TEST_SUITE(ParticleScheduleHierarchyTest, ParticlesDatas);


TYPED_TEST(ParticleScheduleHierarchyTest, testing_inject_ghost_layer)
{
    using ParticleArray_t             = TypeParam::TestParam::PhareTypes::ParticleArray_t;
    using GridLayout_t                = TypeParam::GridLayout_t;
    auto constexpr static dim         = TypeParam::dim;
    auto constexpr static ghost_cells = GridLayout_t::nbrParticleGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            pop.domainParticles().clear();
            EXPECT_EQ(pop.domainParticles().size(), 0);
            if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS))
                for (auto const& tile : pop.domainParticles()())
                {
                    EXPECT_EQ(tile().size(), 0);
                }

            GridLayout_t const layout{phare_box_from<dim>(patch->getBox())};
            auto const ghostBox = grow(layout.AMRBox(), ghost_cells);
            for (auto const& box : ghostBox.remove(layout.AMRBox()))
                core::add_ghost_particles(pop.patchGhostParticles(), layout, box, ppc);
        }
        rm.setTime(ions, *patch, 1);
    }

    this->hierarchy.messenger->fillIonGhostParticles(ions, *lvl0, 0);

    auto n_ghost_cells_for_neighbours = [&](auto& patch) {
        auto domainSamBox    = patch->getBox();
        auto const domainBox = phare_box_from<dim>(domainSamBox);
        return core::sum_from(
            core::generate_from(
                [](auto const& el) { return phare_box_from<dim>(el); },
                SAMRAI::hier::HierarchyNeighbors{*this->hierarchy.basicHierarchy->hierarchy(),
                                                 patch->getPatchLevelNumber(),
                                                 patch->getPatchLevelNumber()}
                    .getSameLevelNeighbors(domainSamBox, patch->getPatchLevelNumber())),
            [&](auto& el) { return (*(grow(el, ghost_cells) * domainBox)).size(); });
    };

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch       = rm.setOnPatch(*patch, ions);
        auto const domainBox   = phare_box_from<dim>(patch->getBox());
        auto const check       = [&](auto const& p) { EXPECT_TRUE(isIn(p, domainBox)); };
        auto const check_array = [&](auto const& array) {
            for (auto const& p : array)
                check(p);
        };
        auto const ncells = n_ghost_cells_for_neighbours(patch);
        for (auto& pop : ions)
        {
            EXPECT_EQ(pop.domainParticles().size(), ncells * ppc);

            if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS))
                for (auto const& tile : pop.domainParticles()())
                    check_array(tile());
            else
                check_array(pop.domainParticles());
        }
    }
}

} // namespace PHARE::amr



int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
