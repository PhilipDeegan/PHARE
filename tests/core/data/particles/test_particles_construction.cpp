

#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_comparator.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>


namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{3});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{10});


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
struct ParticleArrayConstructionTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static impl        = Param::impl;

    using GridLayout_t    = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;
    using ParticleArray_t = ParticleArray<
        dim, ParticleArrayInternals<dim, layout_mode, StorageMode::VECTOR, alloc_mode, impl>>;

    GridLayout_t layout{cells};
};



template<typename ParticleArrayConstructionTest_t>
auto run(ParticleArrayConstructionTest_t& self)
{
    // abort_if(self.periodic_neighbours_for(13).size());
}

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>
    // ,TestParam<3, LayoutMode::SoA, AllocatorMode::CPU>
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>

// PHARE_WITH_THRUST(
//     ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::CPU>
// )

PHARE_WITH_GPU(
    ,TestParam<3, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED, 2>
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, 2>
)

>;
// clang-format on



TYPED_TEST_SUITE(ParticleArrayConstructionTest, Permutations_t, );

TYPED_TEST(ParticleArrayConstructionTest, test_move_swap_etc)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto levelGhostParticles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticles, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesNew = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesNew, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesOld = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesOld, this->layout.AMRBox(), ppc);

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    levelGhostParticlesNew.clear();
    levelGhostParticles.clear();
    levelGhostParticles = levelGhostParticlesOld;

    EXPECT_EQ(levelGhostParticlesNew.size(), 0);
    EXPECT_EQ(levelGhostParticlesOld.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(levelGhostParticles.size(), ppc * std::pow(cells, TestFixture::dim));
}

TYPED_TEST(ParticleArrayConstructionTest, test_copy_move)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(particles, this->layout.AMRBox(), ppc);

    auto copy       = particles;
    auto const move = std::move(particles);

    EXPECT_EQ(move.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(copy.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_TRUE(compare_particles(copy, move));
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
