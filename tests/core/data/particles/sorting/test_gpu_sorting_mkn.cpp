

#include "core/vector.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"

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
        ps_cpu.push_back(Particle_t{0.01, 1., iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
        ps_gpu.push_back(Particle_t{0.01, 1., iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
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
        ParticleArraySorter<ParticleArray_cpu>{ps_cpu, domain_box}();
        ParticleArraySorter<ParticleArray_gpu>{ps_gpu, domain_box}();
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


using AoSGPUParticleArray
    = ParticleArray<3, ParticleArrayInternals<3, LayoutMode::AoS, StorageMode::VECTOR,
                                              PHARE::AllocatorMode::GPU_UNIFIED>>;

using SoAGPUParticleArray
    = ParticleArray<3, ParticleArrayInternals<3, LayoutMode::SoA, StorageMode::VECTOR,
                                              PHARE::AllocatorMode::GPU_UNIFIED>>;


using ParticlesArrays
    = testing::Types<ParticleArrayFixture3D<AoSParticleArray<3>, AoSGPUParticleArray>/*,
                     ParticleArrayFixture3D<SoAParticleArray<3>, SoAGPUParticleArray>*/ //
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


struct K
{
    K(int argc, char** argv) { Kokkos::initialize(argc, argv); }
    ~K() { Kokkos::finalize(); }
};


int main(int argc, char** argv)
{
    K k{argc, argv};

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
