
#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_sorter.hpp"

#include "tests/core/data/particles/test_particles.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

template<typename ParticleArray_cpu, typename ParticleArray_gpu>
class ParticleArrayFixture3D
{
    static constexpr std::size_t ndim = 3;
    using Particle_t                  = Particle<ndim>;
    using box_t                       = Box<int, ndim>;

public:
    void add(std::array<int, ndim> const& iCell)
    {
        ps_cpu.push_back(Particle_t{0.01, iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
        ps_gpu.push_back(Particle_t{0.01, iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
    }

    void print()
    {
        for (std::size_t i = 0; i < ps_cpu.size(); ++i)
            std::cout << __LINE__ << " " << Point{ps_cpu.iCell(i)} << " "
                      << flattener(ps_cpu.iCell(i)) << std::endl;
        std::cout << std::endl;


        for (std::size_t i = 0; i < ps_gpu.size(); ++i)
            std::cout << __LINE__ << " " << Point{ps_gpu.iCell(i)} << " "
                      << flattener(ps_gpu.iCell(i)) << std::endl;
        std::cout << std::endl;
    }

    void sort()
    {
        // mkn::gpu::Pointer c{ps_cpu.data()};
        // mkn::gpu::Pointer g{ps_gpu.data()};

        // PHARE_LOG_LINE_SS(c.is_host_ptr());

        // PHARE_LOG_LINE_SS(g.is_host_ptr());
        // PHARE_LOG_LINE_SS(g.is_managed_ptr());

        sort_particles(ps_cpu, domain_box);
        sort_particles(ps_gpu, domain_box);
        // ParticleArraySorter<ParticleArray_cpu>{*ps_cpu, domain_box}();
        // ParticleArraySorter<ParticleArray_gpu>{*ps_gpu, domain_box}();
    }

    ParticleArray_cpu ps_cpu;
    ParticleArray_gpu ps_gpu;
    box_t domain_box{{0, 0, 0}, {2, 2, 2}};
    CellFlattener<box_t> flattener{domain_box};
};

template<typename ParticlesData>
struct ParticleArraySortingTest : public ::testing::Test
{
};



using AoSGPUParticleArray = ParticleArray<ParticleArrayOptions{
    3, LayoutMode::AoS, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>;


using SoAGPUParticleArray = ParticleArray<ParticleArrayOptions{
    3, LayoutMode::SoA, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>;


// SoA: by_delta walks the array flat, which SoA no longer supports
using ParticlesArrays
    = testing::Types<ParticleArrayFixture3D<AoSParticleArray<3>, AoSGPUParticleArray> //
                     >;

TYPED_TEST_SUITE(ParticleArraySortingTest, ParticlesArrays, );

TYPED_TEST(ParticleArraySortingTest, _3d_sorting_test)
{
    std::size_t constexpr static N = 5;

    TypeParam fixture{};

    std::vector<std::array<int, 3>> const iCells{{2, 2, 2}, {2, 1, 0}, {1, 1, 1}, //
                                                 {1, 0, 1}, {0, 2, 0}, {0, 0, 0},
                                                 {1, 1, 2}, {0, 1, 0}, {2, 0, 2}};
    for (std::size_t j = 0; j < N; ++j)
        for (auto const& iCell : iCells)
            fixture.add(iCell);

    ASSERT_EQ(fixture.ps_cpu.size(), fixture.ps_gpu.size());

    fixture.sort();

    std::vector<std::array<int, 3>> const expected{{0, 0, 0}, {0, 1, 0}, {0, 2, 0}, //
                                                   {1, 0, 1}, {1, 1, 1}, {1, 1, 2},
                                                   {2, 0, 2}, {2, 1, 0}, {2, 2, 2}};

    ASSERT_EQ(fixture.ps_cpu.size(), expected.size() * N);

    for (std::size_t e = 0, i = 0; i < fixture.ps_cpu.size(); ++e, i += N)
        for (std::size_t j = 0; j < N; ++j)
            ASSERT_EQ(fixture.ps_cpu.iCell(i + j), expected[e]);

    for (std::size_t e = 0, i = 0; i < fixture.ps_gpu.size(); ++e, i += N)
        for (std::size_t j = 0; j < N; ++j)
            ASSERT_EQ(fixture.ps_gpu.iCell(i + j), expected[e]);
}


// tiled AoS: each tile sorted on its own by (cell, delta), the GPU sort must match the CPU one
TEST(ParticleArraySortingTest, _3d_AoSTS_sorting_test)
{
    std::size_t constexpr static dim = 3;
    using CPU_t = ParticleArray<
        ParticleArrayOptions{dim, LayoutMode::AoSTS, StorageMode::VECTOR, PHARE::AllocatorMode::CPU}>;
    using GPU_t = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoSTS, StorageMode::VECTOR,
                                                     PHARE::AllocatorMode::GPU_UNIFIED}>;

    Box<int, dim> const box{{0, 0, 0}, {15, 15, 15}}; // several tiles per dimension
    std::size_t constexpr ghost_cells = 2;
    auto cpu = make_particles<CPU_t>(box, ghost_cells);
    auto gpu = make_particles<GPU_t>(box, ghost_cells);
    add_particles(cpu, box, 5);
    add_particles(gpu, box, 5);
    ASSERT_GT(cpu().size(), 1ul);

    // same random deltas and the same shuffled order per tile on both
    for (std::size_t t = 0; t < cpu().size(); ++t)
    {
        auto& c = cpu()[t]();
        auto& g = gpu()[t]();
        delta_disperse(c, t);
        delta_disperse(g, t);
        std::shuffle(c.begin(), c.end(), std::mt19937_64{t});
        std::shuffle(g.begin(), g.end(), std::mt19937_64{t});
    }

    sort_particles(cpu, box);
    sort_particles(gpu, box);

    CellFlattener<Box<int, dim>> const cf{box};
    for (std::size_t t = 0; t < cpu().size(); ++t)
    {
        auto const& c = cpu()[t]();
        auto const& g = gpu()[t]();
        ASSERT_EQ(c.size(), g.size());
        for (std::size_t i = 0; i < c.size(); ++i)
        {
            EXPECT_EQ(c[i].iCell(), g[i].iCell());
            EXPECT_EQ(c[i].delta(), g[i].delta());
            if (i > 0)
            {
                auto const prev = cf(c[i - 1].iCell()), here = cf(c[i].iCell());
                ASSERT_LE(prev, here);
                if (prev == here)
                    ASSERT_LE(as_tuple(c[i - 1].delta()), as_tuple(c[i].delta()));
            }
        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
