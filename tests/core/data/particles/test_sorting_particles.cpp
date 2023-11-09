
#include <random>

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_sorting.hpp"
#include "core/data/particles/particle_utilities.hpp"
#include "core/utilities/box/box.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

std::size_t constexpr static nbrParticles = 55;

template<typename ParticleArray_>
struct SortingTest : public ::testing::Test
{
    using ParticleArray_t                          = ParticleArray_;
    using box_t                                    = typename ParticleArray_t::box_t;
    auto constexpr static dim                      = ParticleArray_t::dimension;
    std::int32_t constexpr static particle_ghosts  = 1;
    std::int32_t constexpr static ghost_safe_layer = 2;

    SortingTest()
        : domainBox{ConstArray<int, dim>(50), ConstArray<int, dim>(54)}
        , particles(grow(domainBox, ghost_safe_layer), nbrParticles) // includes ghosts
    {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> dist(particles.box().lower[0], particles.box().upper[0]);

        for (std::size_t iPart = 0; iPart < nbrParticles; ++iPart)
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
                particles[iPart].iCell[iDim] = dist(gen);
    }

    box_t domainBox;
    box_t ghostBox{grow(domainBox, particle_ghosts)};
    ParticleArray_t particles;
};

using ParticleList = testing::Types<ParticleArray<1>, ParticleArray<2> /*, ParticleArray<3>*/>;

TYPED_TEST_SUITE(SortingTest, ParticleList);

TYPED_TEST(SortingTest, TestRangeFinder)
{
    auto static constexpr dim = TestFixture::dim;
    using box_t               = typename TestFixture::box_t;

    CountingSort<TypeParam, TypeParam::dimension> counting_sort;

    counting_sort.setup(nbrParticles, this->particles.box());

    ParticleCountSorting<TypeParam> sorter{this->particles, counting_sort};
    sorter(); //

    box_t find_box{ConstArray<int, dim>(52), ConstArray<int, dim>(53)};
    ParticleCountRangeFinder<TypeParam> finder{this->particles};

    std::uint32_t count_ppc = 0;
    for (auto const& cell : this->particles.box())
        if (isIn(cell, find_box))
            count_ppc += this->particles.ppc(cell);
    EXPECT_EQ(count_ppc,
              sum_from([](auto const& range) { return range.size(); }, finder.ranges(find_box)));

    count_ppc = 0;
    for (auto const& cell : this->particles.box())
        count_ppc += this->particles.ppc(cell);

    EXPECT_EQ(this->particles.size(), count_ppc);
    EXPECT_EQ(this->particles.size(), sum_from([](auto const& range) { return range.size(); },
                                               finder.ranges(this->particles.box())));

    EXPECT_EQ(this->particles, finder.copy(this->particles.box()));
}

TYPED_TEST(SortingTest, TestRangeFinderSelectGhostLayer)
{
    auto static constexpr dim = TestFixture::dim;
    using ParticleFinder      = ParticleCountRangeFinder<TypeParam>;

    CountingSort<TypeParam, TypeParam::dimension> counting_sort;
    counting_sort.setup(nbrParticles, this->particles.box());

    ParticleCountSorting<TypeParam> sorter{this->particles, counting_sort};
    sorter(); // sort / count ppc / map ppc via icell
    ParticleFinder finder{this->particles};

    for (auto outside_ghost_box : this->particles.box().remove(this->ghostBox))
        for (auto const& range : finder.ranges(outside_ghost_box))
            for (auto const& p : range)
                EXPECT_TRUE(not isIn(p.iCell, this->ghostBox));

    ParticleRangeEraser<TypeParam> eraser{sorter};
    eraser.erase_outside(this->ghostBox);

    for (auto const& p : sorter.particles)
        EXPECT_EQ(isIn(p.iCell, this->ghostBox), true);


    sorter(); // to select domain after erase
    for (auto& range : finder.ranges(this->domainBox))
        for (auto const& p : range)
            EXPECT_TRUE(isIn(p.iCell, this->domainBox));

    eraser.erase(this->domainBox);

    for (auto const& p : this->particles)
    {
        bool inDomainBox = isIn(p.iCell, this->domainBox);
        bool inGhostBox  = isIn(p.iCell, this->ghostBox);
        EXPECT_EQ(true, not inDomainBox and inGhostBox);
    }

    sorter(); // to select ghost layer only
    for (auto const& range : finder.ranges(this->particles.box()))
    {
        for (auto const& p : range)
        {
            bool inDomainBox = isIn(p.iCell, this->domainBox);
            bool inGhostBox  = isIn(p.iCell, this->ghostBox);
            EXPECT_EQ(true, not inDomainBox and inGhostBox);
        }
    }
}

TYPED_TEST(SortingTest, TestRangeFinderSelectDomain)
{
    auto static constexpr dim = TestFixture::dim;
    using ParticleFinder      = ParticleCountRangeFinder<TypeParam>;

    CountingSort<TypeParam, TypeParam::dimension> counting_sort;
    counting_sort.setup(nbrParticles, this->particles.box());

    ParticleCountSorting<TypeParam> sorter{this->particles, counting_sort};
    sorter(); // sort / count ppc / map ppc via icell
    ParticleFinder finder{this->particles};

    for (auto outside_domain_box : this->particles.box().remove(this->domainBox))
        for (auto const& range : finder.ranges(outside_domain_box))
            for (auto const& p : range)
                EXPECT_TRUE(not isIn(p.iCell, this->domainBox));

    ParticleRangeEraser<TypeParam>{sorter}.erase_outside(this->domainBox);

    for (auto const& p : sorter.particles)
        EXPECT_EQ(isIn(p.iCell, this->domainBox), true);

    sorter(); // to select domain after erase

    for (auto& range : finder.ranges(this->domainBox))
        for (auto const& p : range)
            EXPECT_TRUE(isIn(p.iCell, this->domainBox));
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
