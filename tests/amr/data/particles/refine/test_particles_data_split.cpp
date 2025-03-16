
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include "amr/data/particles/refine/particles_data_split.hpp"
#include "amr/data/particles/particles_data.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "tests/amr/amr.hpp"

#include "gtest/gtest.h"
#include <SAMRAI/pdat/CellGeometry.h>

using namespace PHARE;
using namespace PHARE::core;
using namespace PHARE::amr;


std::size_t constexpr static dim            = 2;
std::size_t constexpr static interp         = 1;
std::size_t constexpr static ghost_cells    = 1;
std::size_t constexpr static nb_split_parts = 4;
std::size_t constexpr static ppc            = 10;
std::size_t constexpr static cells          = 4;
using GridLayout_t = TestGridLayout<typename core::PHARE_Types<dim, interp>::GridLayout_t>;
using Box_t        = PHARE::core::Box<int, dim>;


auto static particles_dict()
{
    initializer::PHAREDict dict;
    dict["tile_size"]    = std::size_t{cells};
    dict["interp_order"] = interp;
    return dict;
}


template<std::size_t _dim, auto lm, auto am, std::uint8_t _impl = 0>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
    auto constexpr static impl        = _impl;
};


template<typename Param>
struct AParticlesData
{
    static constexpr auto dim         = Param::dim;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    using ParticleArray_t = ParticleArray<
        dim, ParticleArrayInternals<dim, layout_mode, StorageMode::VECTOR, alloc_mode, /*impl*/ 2>>;
};


template<typename AParticlesData>
struct Patch
{
    SAMRAI::tbox::Dimension static inline const dimension{dim};
    SAMRAI::hier::BlockId static inline const blockId{0};
    SAMRAI::hier::IntVector static inline const ghostVec{dimension, ghost_cells};

    using ParticleArray_t = AParticlesData::ParticleArray_t;
    using ParticlesData_t = ParticlesData<ParticleArray_t>;


    Patch(GridLayout_t const& _layout)
        : layout{_layout}
    {
    }


    GridLayout_t const layout;
    SAMRAI::hier::Box const domain{samrai_box_from(layout.AMRBox())};
    std::shared_ptr<SAMRAI::hier::BoxGeometry> const geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(domain, ghostVec)};

    std::shared_ptr<ParticlesData<ParticleArray_t>> data{
        std::make_shared<ParticlesData<ParticleArray_t>>(domain, ghostVec, "name", [&]() {
            return make_particles<ParticleArray_t>(*layout, particles_dict());
        })};

    SAMRAI::hier::Box const mask{data->getGhostBox()};
};


template<typename ParticlesData>
struct ParticlesDataTest : public ::testing::Test, public ParticlesData
{
    using ParticleArray_t = ParticlesData::ParticleArray_t;

    using Splitter
        = PHARE::amr::Splitter<PHARE::core::DimConst<dim>, PHARE::core::InterpConst<interp>,
                               PHARE::core::RefinedParticlesConst<nb_split_parts>>;

    using Refiner = ParticlesRefining<ParticleArray_t, ParticlesDataSplitType::interior, Splitter>;

    ParticlesDataTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
            {
                auto const cellx = i * cells;
                auto const celly = j * cells;
                patches.emplace_back(Box_t{Point{cellx, celly}, Point{cellx + off, celly + off}});
            }

        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));

        for (auto& patch : patches)
            add_particles_in(patch.data->domainParticles, patch.layout.AMRBox(), ppc);
    }


    std::vector<Patch<ParticlesData>> patches;
    Patch<ParticlesData> L1{GridLayout_t{Box_t{Point{2, 4}, Point{17, 19}}}};
};


// clang-format off
using ParticlesDatas = testing::Types< //

    AParticlesData<TestParam<dim, LayoutMode::AoSMapped, AllocatorMode::CPU>>
   ,AParticlesData<TestParam<dim, LayoutMode::AoSTS, AllocatorMode::CPU>>
   // ,AParticlesData<TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>>

PHARE_WITH_GPU(

   // ,AParticlesData<TestParam<3, LayoutMode::SoA, AllocatorMode::CPU>>
   // ,AParticlesData<TestParam<3, LayoutMode::SoATS, AllocatorMode::CPU>>

)

>;
// clang-format on

TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas);

namespace PHARE::amr
{


TYPED_TEST(ParticlesDataTest, splitWorks)
{
    auto& dst = this->L1;

    for (auto const& src : this->patches)
        typename TestFixture::Refiner{*src.data, *dst.data}.forBoxes(
            std::array{dst.layout.AMRBox()});

    std::size_t const expected = this->L1.layout.AMRBox().size() * ppc;
    EXPECT_EQ(expected, dst.data->domainParticles.size());

    // if constexpr (TestFixture::ParticleArray_t::layout_mode == LayoutMode::AoSTS)
    // {
    //     for (auto const& tile : dst.data->domainParticles())
    //     {
    //         for (auto const& particle : tile())
    //         {
    //             // EXPECT_TRUE(isIn(particle, tile));
    //             PHARE_LOG_LINE_SS(Point{particle.iCell()} << " " << *tile << " "
    //                                                       << isIn(particle, tile));
    //         }
    //     }
    // }
    // else
    // {
    //     for (auto const& particle : dst.data->domainParticles)
    //     {
    //         EXPECT_TRUE(isIn(particle, dst.layout.AMRBox()));
    //     }
    // }
}


} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
