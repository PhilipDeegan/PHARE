
#include "phare_core.hpp"

#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>
#include "core/data/grid/grid_tiles.hpp"
#include <core/data/ndarray/ndarray_vector.hpp>

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"


#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>


#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 100;


template<SimOpts opts>
struct TestParam
{
    auto constexpr static dim    = opts.dimension;
    auto constexpr static interp = opts.interp_order;

    using PhareTypes       = core::PHARE_Types<opts>;
    using Field_t          = PhareTypes::Field_t;
    using GridLayout_t     = PhareTypes::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
    using Box_t            = core::Box<int, dim>;
    using ParticleArray_t  = PhareTypes::ParticleArray_t;

    using Hierarchy_t = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct FieldScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    using Hierarchy_t         = typename TestParam::Hierarchy_t;
    using ResourceManager_t   = typename Hierarchy_t::ResourcesManagerT;
    auto constexpr static dim = TestParam::dim;

    std::string const configFile
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_config.txt";
    Hierarchy_t hierarchy{configFile};
};


// clang-format off
using FieldDatas = testing::Types<
    TestParam<SimOpts{}>
   ,TestParam<SimOpts{2}>
PHARE_WITH_MKN_GPU(
   ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
   ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
)

>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas, );


TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_schedules)
{
    auto constexpr static dim = TypeParam::dim;
    using TestParam           = TestFixture::TestParam;
    using ParticleArray_t     = TestParam::ParticleArray_t;
    using GridLayout_t        = TestParam::GridLayout_t;

    using FieldData_t = TestFixture::ResourceManager_t::UserField_t::patch_data_type;
    using Interpolating_t
        = core::Interpolating<ParticleArray_t, TestParam::interp, /*atomic_interp*/ false>;


    auto constexpr static interp      = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    auto proton_flux_x_id = *rm.getID("protons_flux_x");


    Interpolating_t interpolate;

    for (auto& patch : *lvl0)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, ions);
        resetMoments(ions);

        core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});

        auto const field_data = SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            patch->getPatchData(proton_flux_x_id));

        for (auto& pop : ions)
            if (pop.name() == "protons")
            {
                PHARE_LOG_LINE_SS(core::sum_field(reduce(pop.density())));
            }
    }

    this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);
    this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);


    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);

        auto const field_data = SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            patch->getPatchData(proton_flux_x_id));

        if (mpi::rank() == 0)
            for (auto& pop : ions)
                if (pop.name() == "protons")
                {
                    auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
                    auto const domGhostBox
                        = layout.AMRGhostBoxFor(pop.flux()[0].physicalQuantity());

                    PHARE_LOG_LINE_SS(core::sum_field(reduce_single(pop.density())));
                }
    }
}




TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_field_refine_schedules)
{
    auto constexpr static dim    = TypeParam::dim;
    using GridLayout_t           = TestFixture::TestParam::GridLayout_t;
    using Field_t                = TestFixture::TestParam::Field_t;
    using FieldData_t            = TestFixture::ResourceManager_t::UserField_t::patch_data_type;
    auto constexpr static interp = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0        = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto lvl1        = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(1);
    auto& rm         = *this->hierarchy.resourcesManagerHybrid;
    auto& electromag = this->hierarchy.hybridModel->state.electromag;
    auto& ions       = this->hierarchy.hybridModel->state.ions;


    this->hierarchy.messenger->fillElectricGhosts(electromag.E, 1, 0);


    auto Bx_id = *rm.getID("EM_B_x");


    for (auto& patch : *lvl1)
    {
        PHARE_LOG_LINE_SS(patch->getGlobalId());
        auto dataOnPatch = rm.setOnPatch(*patch, *this->hierarchy.hybridModel);

        for (auto& pop : ions)
        {
            assert(pop.momentumTensor().isUsable());
        }
        ions.computeFullMomentumTensor();

        auto const field_data = SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            patch->getPatchData(Bx_id));

        if (mpi::rank() == 0)
        {
            auto const Bx = electromag.B[0]; // primal in x/1d

            // PHARE_LOG_LINE_SS(Bx);

            if constexpr (core::is_field_tile_set_v<Field_t>)
            {
                // PHARE_LOG_LINE_SS(reduce_single(field_data->field));
            }
            auto const Ey = electromag.E[1]; // primal in x/1d

            // PHARE_LOG_LINE_SS(Ey);

            if constexpr (core::is_field_tile_set_v<Field_t>)
            {
                // PHARE_LOG_LINE_SS(reduce_single(Ey));
            }
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
