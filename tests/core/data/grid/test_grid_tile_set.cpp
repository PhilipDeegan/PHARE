
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "phare_core.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"


using namespace PHARE;
using namespace PHARE::core;


std::size_t constexpr static interp = 1;
std::size_t constexpr static cells  = 14;


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    auto constexpr static dim = _dim;

    auto constexpr static SimOpts opts{dim, interp, lm, am};

    using PhareTypes = core::PHARE_Types<opts>;

    using Box_t        = PHARE::core::Box<int, dim>;
    using GridLayout_t = TestGridLayout<typename PhareTypes::GridLayout_t>;
};




template<typename TestParam>
struct Patch
{
    auto constexpr static dim = TestParam::dim;

    using PhareTypes   = TestParam::PhareTypes;
    using GridLayout_t = TestParam::GridLayout_t;
    using Box_t        = TestParam::Box_t;

    using Field_t = PhareTypes::Field_t;
    using Grid_t  = PhareTypes::Grid_t;

    Patch(Box_t const& box)
        : layout{box}
    {
    }
    Patch(GridLayout_t const& _layout)
        : layout{_layout}
    {
    }


    GridLayout_t const layout;
    Grid_t rho{"rho", layout, HybridQuantity::Scalar::rho};
    Field_t& rho_v = *rho;
};

template<std::size_t dim, typename TestParam>
struct AGridFieldTest;

template<typename TestParam>
struct AGridFieldTest<1, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
    {
        patches.reserve(3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{Point{cellx}, Point{cellx + off}});
        }
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{3}, Point{12}}};
};


template<typename TestParam>
struct AGridFieldTest<2, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
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
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{3, 3}, Point{12, 12}}};
};


template<typename TestParam>
struct AGridFieldTest<3, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    auto const cellx = i * cells;
                    auto const celly = j * cells;
                    auto const cellz = k * cells;
                    patches.emplace_back(Box_t{Point{cellx, celly, cellz},
                                               Point{cellx + off, celly + off, cellz + off}});
                }
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{3, 3, 3}, Point{12, 12, 12}}};
};


template<typename TestParam>
struct GridFieldTest : public ::testing::Test, public AGridFieldTest<TestParam::dim, TestParam>
{
    using Super        = AGridFieldTest<TestParam::dim, TestParam>;
    using GridLayout_t = TestParam::GridLayout_t;

    using Super::L1;
    using Super::patches;

    GridFieldTest() {}
};


// clang-format off
using ParticlesDatas = testing::Types< //

PHARE_WITH_THRUST(
    TestParam<1, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>

PHARE_WITH_GPU(
   ,TestParam<2, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
)

)

>;
// clang-format on

TYPED_TEST_SUITE(GridFieldTest, ParticlesDatas);


namespace PHARE::core
{

TYPED_TEST(GridFieldTest, test_compiles)
{
    auto& L1 = this->L1;

    EXPECT_EQ(L1.rho().data(), L1.rho_v().data());

    for (auto& tile : L1.rho())
    {
        Point const lower{tile.lower};
        EXPECT_EQ((*L1.rho).at(lower), L1.rho_v.at(lower));
    }
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
