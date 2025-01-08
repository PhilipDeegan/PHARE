#ifndef PHARE_CORE_UTILITIES_KERNELS_DETAIL_DEF_LAUNCHERS_HPP
#define PHARE_CORE_UTILITIES_KERNELS_DETAIL_DEF_LAUNCHERS_HPP

#include "core/def.hpp"

#if PHARE_HAVE_MKN_GPU

namespace PHARE::kernel::gpu
{
struct Params
{
    dim3 block_size, grid_size;
};
} // namespace PHARE::kernel::gpu

namespace PHARE::kernel
{

auto inline warp_size() _PHARE_ALL_FN_ // make runtime interrogation?
{
#if !defined(MKN_GPU_CUDA) || !defined(MKN_GPU_ROCM)
    throw std::runtime_error("warpsize is not available for cpu");
#elif defined(MKN_GPU_ROCM)
    return 32ull;
#elif defined(MKN_GPU_CUDA)
    return 32ull;
#endif
}

auto inline thread_idx() _PHARE_DEV_FN_
{
    if (CompileOptions::WithGpu)
    {
        return PHARE_WITH_GPU(threadIdx.x);
    }
    // else compile error
}

auto inline block_idx() _PHARE_DEV_FN_
{
    if (CompileOptions::WithGpu)
    {
        return PHARE_WITH_GPU(blockIdx.x);
    }
    // else compile error
}

auto inline idx() _PHARE_DEV_FN_ // global
{
    if (CompileOptions::WithMknGpu)
    {
        return PHARE_WITH_MKN_GPU(mkn::gpu::idx());
    }
    // else compile error
}

} // namespace PHARE::kernel


#endif // PHARE_HAVE_MKN_GPU

#endif /* PHARE_CORE_UTILITIES_KERNELS_DETAIL_DEF_LAUNCHERS_HPP */
