// SoAVXChunked: counting sort by cell + per cell offsets, box selection as particle ranges

#include "core/def/phare_config.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"

#include "gtest/gtest.h"

#include <map>
#include <random>
#include <vector>
#include <cstddef>

namespace PHARE::core
{

template<std::size_t _dim, auto am>
struct TestParam
{
    auto constexpr static dim        = _dim;
    auto constexpr static alloc_mode = am;
};

template<typename Param>
struct SoAVXChunkedTest : public ::testing::Test
{
    auto constexpr static dim        = Param::dim;
    auto constexpr static alloc_mode = Param::alloc_mode;

    using Box_t      = Box<int, dim>;
    using Particle_t = Particle<dim>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, LayoutMode::SoAVXChunked, StorageMode::VECTOR,
                                             alloc_mode}>;

    static Box_t make_box(int const lo, int const up)
    {
        Box_t box;
        for (std::size_t d = 0; d < dim; ++d)
        {
            box.lower[d] = lo;
            box.upper[d] = up;
        }
        return box;
    }

    // particles in random cells of the shell (box + 1), weight = unique id
    SoAVXChunkedTest()
        : particles{make_box(0, 9)}
    {
        std::mt19937 gen{1337};
        auto const shell = particles.shell();
        std::uniform_int_distribution<int> cell_dist{shell.lower[0], shell.upper[0]};
        std::uniform_real_distribution<double> real_dist{0, 1};

        std::size_t const n = shell.size() * 10;
        for (std::size_t i = 0; i < n; ++i)
        {
            Particle_t p;
            p.weight() = static_cast<double>(i);
            for (std::size_t d = 0; d < dim; ++d)
            {
                p.iCell()[d] = cell_dist(gen);
                p.delta()[d] = real_dist(gen);
            }
            for (std::size_t d = 0; d < 3; ++d)
                p.v()[d] = real_dist(gen);
            originals.push_back(p);
            particles.push_back(p);
        }
    }

    std::vector<Particle_t> originals;
    ParticleArray_t particles;
};

// clang-format off
using Permutations_t = testing::Types<
    TestParam<1, AllocatorMode::CPU>
   ,TestParam<2, AllocatorMode::CPU>
   ,TestParam<3, AllocatorMode::CPU>
PHARE_WITH_THRUST( // the GPU sort is thrust (not on ROCm yet)
   ,TestParam<1, AllocatorMode::GPU_UNIFIED>
   ,TestParam<2, AllocatorMode::GPU_UNIFIED>
   ,TestParam<3, AllocatorMode::GPU_UNIFIED>
)
>;
// clang-format on

TYPED_TEST_SUITE(SoAVXChunkedTest, Permutations_t, );

TYPED_TEST(SoAVXChunkedTest, sort_keeps_particles_and_builds_offsets)
{
    auto& particles = this->particles;
    EXPECT_FALSE(particles.sorted());

    particles.sort();
    ASSERT_TRUE(particles.sorted());
    ASSERT_EQ(particles.size(), this->originals.size());

    // every particle once, all its components moved with it
    std::vector<bool> seen(particles.size(), false);
    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        auto const p  = particles.copy(i);
        auto const id = static_cast<std::size_t>(p.weight());
        ASSERT_LT(id, seen.size());
        EXPECT_FALSE(seen[id]);
        seen[id] = true;
        EXPECT_EQ(p, this->originals[id]);
    }

    // each cell's range holds exactly its particles, ranges tile [0, size)
    std::size_t expected_start = 0;
    for_each_point(particles.shell(), [&](auto const& cell) {
        auto const [start, size] = particles.cell_range(cell);
        EXPECT_EQ(start, expected_start);
        for (std::size_t i = start; i < start + size; ++i)
        {
            EXPECT_EQ(particles.iCell(i), cell.toArray());
        }
        expected_start += size;
    });
    EXPECT_EQ(expected_start, particles.size());

    particles.push_back(this->originals[0]);
    EXPECT_FALSE(particles.sorted());
}

TYPED_TEST(SoAVXChunkedTest, select_box_ranges)
{
    using Box_t     = TestFixture::Box_t;
    auto& particles = this->particles;
    particles.sort();

    auto const box   = particles.box();
    auto const shell = particles.shell();
    auto const inner = TestFixture::make_box(2, 6);

    // the ghost layer of the shell, as an exporter would take leavers
    auto const layer = shell.remove(box);

    for (auto const& select : std::vector<Box_t>{box, shell, inner})
    {
        std::size_t brute = 0;
        for (std::size_t i = 0; i < particles.size(); ++i)
            brute += isIn(Point{particles.iCell(i)}, select);

        std::size_t ranged = 0, prev_end = 0, n_ranges = 0;
        particles.for_each_range(select, [&](auto const start, auto const size) {
            EXPECT_GE(start, prev_end); // ascending, disjoint
            prev_end = start + size;
            ++n_ranges;
            for (std::size_t i = start; i < start + size; ++i)
            {
                EXPECT_TRUE(isIn(Point{particles.iCell(i)}, select));
            }
            ranged += size;
        });
        EXPECT_EQ(ranged, brute);
        EXPECT_EQ(particles.count(select), brute);
        if (select == shell)
        {
            EXPECT_EQ(n_ranges, 1u); // the whole shell is one range
        }
    }

    std::size_t layer_count = 0;
    for (auto const& b : layer)
        layer_count += particles.count(b);
    EXPECT_EQ(layer_count + particles.count(box), particles.size());

    // a view (span) reads the same offsets
    auto const view = particles.view();
    EXPECT_TRUE(view.sorted());
    EXPECT_EQ(view.count(inner), particles.count(inner));
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
