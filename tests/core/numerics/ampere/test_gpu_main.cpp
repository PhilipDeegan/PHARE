//
//

#include "core/numerics/ampere/ampere.hpp"
#include "core/data/grid/gridlayout_utils.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"



// template<typename Impl>
// struct sycl::is_device_copyable<PHARE::core::GridLayout<Impl>> : std::true_type
// {
// };

// template<std::size_t dim, typename PhysicalQuantity, typename Data_t, auto alloc_mode>
// struct sycl::is_device_copyable<PHARE::core::Field<dim, PhysicalQuantity, Data_t, alloc_mode>>
//     : std::true_type
// {
// };

// template<typename Field_t, typename PhysicalQuantity, std::size_t rank_>
// struct sycl::is_device_copyable<PHARE::core::TensorField<Field_t, PhysicalQuantity, rank_>>
//     : std::true_type
// {
// };



struct AmpereTest : public ::testing::Test
{
};

template<auto alloc_mode, typename GridLayout_t>
auto evolve(GridLayout_t const& layout)
{
    using namespace PHARE::core;
    std::size_t constexpr static dim = GridLayout_t::dimension;
    UsableElectromag<dim, alloc_mode> em{layout};
    UsableVecField<dim, alloc_mode> J{"J", layout, HybridQuantity::Vector::J};
    Ampere_ref<GridLayout_t>{layout}(*em.B, *J);
    return J;
}


template<std::size_t dim, std::size_t cells = 30>
void test()
{
    using PHARE_Types = PHARE::core::PHARE_Types<dim, /*interp=*/3>;
    TestGridLayout<typename PHARE_Types::GridLayout_t> layout{cells};
    EXPECT_EQ(evolve<PHARE::AllocatorMode::CPU>(*layout),
              evolve<PHARE::AllocatorMode::GPU_UNIFIED>(*layout));
}

// TEST(AmpereTest, worksOnGPU_1d)
// {
//     test<1>();
// }

TEST(AmpereTest, worksOnGPU_2d)
{
    test<2>();
}

// TEST(AmpereTest, worksOnGPU_3d)
// {
//     test<3>();
// }

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
