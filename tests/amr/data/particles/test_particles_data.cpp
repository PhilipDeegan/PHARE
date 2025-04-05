
#include "core/utilities/types.hpp"
#include "phare_core.hpp"

// #include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/particles_data.hpp"

// #include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
// #include "tests/core/data/electromag/test_electromag_fixtures.hpp"
// #include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
// #include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"
#include "tests/core/data/particles/test_particles.hpp"


// #include "tests/core/utilities/box/test_box_fixtures.hpp"
#include "tests/amr/amr.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <SAMRAI/pdat/CellGeometry.h>

using namespace PHARE;
using namespace PHARE::core;
using namespace PHARE::amr;


std::size_t constexpr static dim         = 3; // only 3d
std::size_t constexpr static interp      = 1;
std::size_t constexpr static ghost_cells = 1;
using GridLayout_t = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;
using Box_t        = PHARE::core::Box<int, dim>;


auto static const cells   = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const ppc     = get_env_as("PHARE_PPC", std::size_t{1});
auto static const premain = []() { return true; }();

auto static particles_dict()
{
    initializer::PHAREDict dict;
    dict["tile_size"]    = std::size_t{cells};
    dict["interp_order"] = interp;
    return dict;
}

auto static make_shift_for(Box_t const& box)
{
    int const span = cells * 3;
    int const mid  = cells * 3 / 2;
    auto shifts
        = for_N<7, for_N_R_mode::make_array>([&](auto i) { return Point<int, 3>{0, 0, 0}; });
    for_N<3>([&](auto i) {
        int const shift = box.upper[i] < mid ? 1 : -1;
        shifts[i][i]    = span * shift;
    });

    shifts[3] = {shifts[0][0], shifts[1][1], 0};
    shifts[4] = {0, shifts[1][1], shifts[2][2]};
    shifts[5] = {shifts[0][0], 0, shifts[2][2]};
    shifts[6] = {shifts[0][0], shifts[1][1], shifts[2][2]};

    return shifts;
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

    auto static overlap(Patch const& src, Patch& dst)
    {
        return overlap(src, dst, ConstArray<int, dim>());
    }

    auto static overlap(Patch const& src, Patch& dst, auto const shift)
    {
        bool constexpr static overwriteInterior{true};

        SAMRAI::hier::IntVector shiftVec{dimension};
        for (std::size_t i = 0; i < dim; ++i)
            shiftVec[i] = -shift[i];

        SAMRAI::hier::Box const srcMask{src.data->getGhostBox()};
        SAMRAI::hier::Box const fillBox{dst.data->getGhostBox()};
        SAMRAI::hier::Transformation const transformation{shiftVec};
        return std::shared_ptr<SAMRAI::pdat::CellOverlap>{
            std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(dst.geom->calculateOverlap(
                *src.geom, srcMask, fillBox, overwriteInterior, transformation))};
    }


    GridLayout_t const layout;
    SAMRAI::hier::Box const domain{samrai_box_from(layout.AMRBox())};
    std::shared_ptr<SAMRAI::hier::BoxGeometry> const geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(domain, ghostVec)};
    std::unique_ptr<ParticlesData<ParticleArray_t>> data{
        std::make_unique<ParticlesData<ParticleArray_t>>(domain, ghostVec, "name", [&]() {
            return make_particles<ParticleArray_t>(*layout, particles_dict());
        })};
    SAMRAI::hier::Box const mask{data->getGhostBox()};
};


template<typename ParticlesData>
struct ParticlesDataTest : public ::testing::Test, public ParticlesData
{
    using ParticleArray_t = ParticlesData::ParticleArray_t;

    ParticlesDataTest()
    {
        // 27 patches of cells**3 in 3**3 config
        patches.reserve(3 * 3 * 3);
        int const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    int const cellx = i * cells;
                    int const celly = j * cells;
                    int const cellz = k * cells;
                    patches.emplace_back(Box_t{Point{cellx, celly, cellz},
                                               Point{cellx + off, celly + off, cellz + off}});
                }

        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));

        for (auto& patch : patches)
            add_particles_in(patch.data->domainParticles, patch.layout.AMRBox(), ppc);
    }

    auto static overlap(Patch<ParticlesData> const& src, Patch<ParticlesData>& dst)
    {
        return Patch<ParticlesData>::overlap(src, dst);
    }
    auto static overlap(Patch<ParticlesData> const& src, Patch<ParticlesData>& dst,
                        auto const shift)
    {
        return Patch<ParticlesData>::overlap(src, dst, shift);
    }

    std::vector<Patch<ParticlesData>> patches;
};


// clang-format off
using ParticlesDatas = testing::Types< //

    AParticlesData<TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>>
   ,AParticlesData<TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>>
   ,AParticlesData<TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>>

PHARE_WITH_GPU(

   ,AParticlesData<TestParam<3, LayoutMode::SoA, AllocatorMode::CPU>>
   ,AParticlesData<TestParam<3, LayoutMode::SoATS, AllocatorMode::CPU>>

)

>;
// clang-format on


namespace PHARE::core::detail::strings
{
constexpr static std::string_view copyWorks                   = "copyWorks,";
constexpr static std::string_view copyAtPeriodicBoundaryWorks = "copyAtPeriodicBoundaryWorks,";
constexpr static std::string_view packWorks                   = "packWorks,";
constexpr static std::string_view packAtPeriodicBoundaryWorks = "packAtPeriodicBoundaryWorks,";
constexpr static std::string_view cma                         = ",";
} // namespace PHARE::core::detail::strings


TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas);

namespace PHARE::core
{
TYPED_TEST(ParticlesDataTest, copyWorks)
{
    auto constexpr function_id
        = join_string_views_v<detail::strings::copyWorks, TestFixture::ParticleArray_t::type_id,
                              detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    std::size_t pid = 13;
    auto& dst       = this->patches[pid];

    for (std::size_t i = 0; i < pid; ++i)
        dst.data->copy(*this->patches[i].data);
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        dst.data->copy(*this->patches[i].data);

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, copyAtPeriodicBoundaryWorks)
{
    auto constexpr function_id
        = join_string_views_v<detail::strings::copyAtPeriodicBoundaryWorks,
                              TestFixture::ParticleArray_t::type_id, detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    std::size_t const pid   = 0;
    auto& dst               = this->patches[pid];
    auto const dst_ghostbox = grow(dst.layout.AMRBox(), 1);

    // non-periodic neighbours
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        if (auto const overlap = dst_ghostbox * this->patches[i].layout.AMRBox())
            dst.data->copy(*this->patches[i].data);

    // periodic neighbours
    for (auto const& shifter : make_shift_for(dst.layout.AMRBox()))
    {
        auto const shift_box = shift(dst.layout.AMRBox(), shifter);
        auto const ghostbox  = grow(shift_box, 1);

        for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
            if (auto const overlap = ghostbox * this->patches[i].layout.AMRBox())
            {
                auto const celloverlap = TestFixture::overlap(this->patches[i], dst, shifter);
                dst.data->copy(*this->patches[i].data, *celloverlap);
            }
    }

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packWorks)
{
    auto constexpr function_id
        = join_string_views_v<detail::strings::packWorks, TestFixture::ParticleArray_t::type_id,
                              detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    std::size_t pid = 13;
    auto& dst       = this->patches[pid];

    auto for_neighbour = [&](auto i) {
        auto const cellOverlap = TestFixture::overlap(this->patches[i], dst);
        auto const& src        = this->patches[i];
        SAMRAI::tbox::MessageStream particlesWriteStream;
        src.data->packStream(particlesWriteStream, *cellOverlap);
        SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        particlesWriteStream.getBufferStart()};

        dst.data->unpackStream(particlesReadStream, *cellOverlap);
    };

    for (std::size_t i = 0; i < pid; ++i)
        for_neighbour(i);

    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        for_neighbour(i);

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packAtPeriodicBoundaryWorks)
{
    auto constexpr function_id
        = join_string_views_v<detail::strings::packAtPeriodicBoundaryWorks,
                              TestFixture::ParticleArray_t::type_id, detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    std::size_t const pid   = 0;
    auto& dst               = this->patches[pid];
    auto const dst_ghostbox = grow(dst.layout.AMRBox(), 1);

    auto for_neighbour = [&](auto i, auto shift) {
        auto const cellOverlap = TestFixture::overlap(this->patches[i], dst, shift);
        auto const& src        = this->patches[i];
        SAMRAI::tbox::MessageStream particlesWriteStream;
        src.data->packStream(particlesWriteStream, *cellOverlap);
        SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        particlesWriteStream.getBufferStart()};

        dst.data->unpackStream(particlesReadStream, *cellOverlap);
    };

    // non-periodic neighbours
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        if (auto const overlap = dst_ghostbox * this->patches[i].layout.AMRBox())
            for_neighbour(i, ConstArray<int, dim>());

    // periodic neighbours
    for (auto const& shifter : make_shift_for(dst.layout.AMRBox()))
    {
        auto const shift_box = shift(dst.layout.AMRBox(), shifter);
        auto const ghostbox  = grow(shift_box, 1);

        for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
            if (auto const overlap = ghostbox * this->patches[i].layout.AMRBox())
                for_neighbour(i, shifter);
    }

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
