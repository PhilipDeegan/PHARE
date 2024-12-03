#ifndef PHARE_CORE_UTILITIES_KERNELS_HPP
#define PHARE_CORE_UTILITIES_KERNELS_HPP


#include "core/def/phare_config.hpp"

#if PHARE_HAVE_GPU

#include "core/utilities/box/box.hpp"
#include "core/utilities/kernels/launchers.hpp"

namespace PHARE::kernel::gpu
{
struct Params
{
    dim3 block_size, grid_size;
};
} // namespace PHARE::kernel::gpu

namespace PHARE::kernel
{

auto inline thread_idx()
{
    if (CompileOptions::WithGpu)
    {
        return PHARE_WITH_GPU(threadIdx.x);
    }
    else
        throw std::runtime_error("launch::thread_idx NO ALTERNATIVE");
}

auto inline block_idx()
{
    if (CompileOptions::WithGpu)
    {
        return PHARE_WITH_GPU(blockIdx.x);
    }
    else
        throw std::runtime_error("launch::block_idx NO ALTERNATIVE");
}

auto inline idx() // global
{
    if (CompileOptions::WithMknGpu)
    {
        return PHARE_WITH_MKN_GPU(mkn::gpu::idx());
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}

template<std::size_t dim, typename... Args>
void launch(core::Box<std::uint32_t, dim> const& box, std::size_t const& size, Args&&... args)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        throw std::runtime_error("launch::kernel NO IMPL");
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(
            { core::gpu::BoxCellNLauncher<core::Box<std::uint32_t, dim>>{box, size}(args...); });
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}

template<typename... Args>
void launch(std::size_t const& size, Args&&... args)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        throw std::runtime_error("launch::kernel NO IMPL");
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU({ mkn::gpu::GDLauncher{size}(args...); });
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}


} // namespace PHARE::kernel

#endif // PHARE_HAVE_GPU

#endif /* PHARE_CORE_UTILITIES_KERNELS_HPP */
