#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode = PHARE::AllocatorMode::CPU,
         std::uint8_t _impl = 0>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), PHARE::AllocatorMode>);
    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
    auto constexpr static impl        = _impl;

    using ParticleArray_t = ParticleArray<
        dim, ParticleArrayInternals<dim, layout_mode, StorageMode::VECTOR, alloc_mode, impl>>;
};

template<typename Param>
class ParticleArrayFixture3D
{
    static constexpr std::size_t ndim = 3;
    using Particle_t                  = Particle<ndim>;
    using box_t                       = Box<int, ndim>;
    using ParticleArray_t             = Param::ParticleArray_t;

public:
    void add(std::array<int, ndim> const& iCell)
    {
        particles.push_back(Particle_t{0.01, 1., iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
    }


    void sort()
    {
        auto view    = particles.view();
        using view_t = decltype(view);
        ParticleArraySorter<view_t>{view, domain_box}();
    }

    void print()
    {
        for (auto const& particle : particles)
            std::cout << __LINE__ << " " << Point{particle.iCell()} << std::endl;
        std::cout << std::endl;
    }

    ParticleArray_t particles;
    box_t domain_box{{0, 0, 0}, {2, 2, 2}};
};

template<typename ParticlesData>
struct ParticleArraySortingTest : public ::testing::Test
{
};



using ParticlesArrays = testing::Types<ParticleArrayFixture3D<TestParam<3, LayoutMode::SoAVX>>,
                                       ParticleArrayFixture3D<TestParam<3, LayoutMode::SoA>>,
                                       ParticleArrayFixture3D<TestParam<3, LayoutMode::AoS>> //
                                       >;

TYPED_TEST_SUITE(ParticleArraySortingTest, ParticlesArrays, );

TYPED_TEST(ParticleArraySortingTest, _3d_sorting_test)
{
    TypeParam fixture{};
    auto& particles = fixture.particles;
    std::vector<std::array<int, 3>> iCells{{2, 2, 2}, {2, 1, 0}, {1, 1, 1}, //
                                           {1, 0, 1}, {0, 2, 0}, {0, 0, 0},
                                           {1, 1, 2}, {0, 1, 0}, {2, 0, 2}};
    for (auto const& iCell : iCells)
        fixture.add(iCell);

    fixture.sort();

    std::vector<std::array<int, 3>> expected{{0, 0, 0}, {0, 1, 0}, {0, 2, 0}, //
                                             {1, 0, 1}, {1, 1, 1}, {1, 1, 2},
                                             {2, 0, 2}, {2, 1, 0}, {2, 2, 2}};
    ASSERT_EQ(particles.size(), expected.size());
    for (std::size_t i = 0; i < particles.size(); ++i)
        ASSERT_EQ(particles.iCell(i), expected[i]);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
