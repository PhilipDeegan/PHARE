// expects to be included
#include <thrust/partition.h>
#include <thrust/execution_policy.h>

// TYPED_TEST(DimConstTest, particle_array_is_host_to_device_copyable)
// {
//     constexpr auto dim = TypeParam{};
//     using Particle_t =  PHARE::core::Particle<dim>
//     using host_t = PHARE::core::ParticleArray<TypeParam{}, typename std::vector<Particle_t>::allocator_type>;
//     using dev_t = PHARE::Vector<Particle_t>;
    
//     host_t host_particles(
// }


TYPED_TEST(VectorTest, is_host_to_device_copyable)
{
    using T = typename TypeParam::value_type;

    std::vector<T> vec0(10, 1);
    auto vec1 = TypeParam::make(10);
    EXPECT_EQ(PHARE::core::sum(vec1), 0);

    TypeParam::copy(vec1, vec0);
    EXPECT_EQ(PHARE::core::sum(vec1), 10);
}

TYPED_TEST(VectorTest, is_device_to_host_copyable)
{
    using T = typename TypeParam::value_type;

    auto vec0 = TypeParam::make(10);
    TypeParam::fill(vec0, 1);
    EXPECT_EQ(PHARE::core::sum(vec0), 10);

    std::vector<T> vec1(10, 0);
    EXPECT_EQ(PHARE::core::sum(vec1), 0);
    TypeParam::copy(vec1, vec0);
    EXPECT_EQ(PHARE::core::sum(vec1), 10);
}

struct is_even
{
    __device__ bool operator()(const int& x) { return (x % 2) == 0; }
};

TYPED_TEST(VectorTest, thrust_device_partition)
{
    auto constexpr N = 1024;
    auto A           = TypeParam::make(N);

    assert(mkn::gpu::Pointer{A.data()}.is_managed_ptr());
    for (std::size_t i = 0; i < N; ++i)
        A[i] = i;

    // there might be other versions that work
    thrust::stable_partition(A.data(), A.data() + N, is_even());

    for (std::size_t i = 0; i < N / 2; ++i)
        EXPECT_EQ(static_cast<int>(A[i]) % 2, 0);
    for (std::size_t i = N / 2; i < N; ++i)
        EXPECT_EQ(static_cast<int>(A[i]) % 2, 1);
}
