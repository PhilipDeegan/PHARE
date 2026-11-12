#ifndef PHARE_CORE_UTILITIES_MEMORY_HPP
#define PHARE_CORE_UTILITIES_MEMORY_HPP


#include "core/def/kokkos_tools.hpp"
#include "core/def/phare_config.hpp"
#include "core/utilities/memory/kokkos_allocator.hpp"
#include "core/utilities/memory/umpire_allocator.hpp"



namespace PHARE::core::gpu
{


template<typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        PHARE_WITH_RAJA(PHARE::core::raja::copy(dst, src, size));
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU({
            assert(
                not(mkn::gpu::Pointer{dst}.is_host_ptr() and mkn::gpu::Pointer{src}.is_host_ptr()));
            mkn::gpu::copy(dst, src, size);
        })
    }
    else
        throw std::runtime_error("Vector::copy NO ALTERNATIVE");
}


template<typename Vector, typename Value>
void fill(Vector& dst, Value const val)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        PHARE_WITH_RAJA(PHARE::core::raja::set(dst.data(), val, dst.size()));
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(mkn::gpu::fill(dst, val));
    }
    else if (CompileOptions::WithKokkos)
    {
        PHARE_WITH_KOKKOS({ //
            auto data = dst.data();
            Kokkos::parallel_for("fill", 5, KOKKOS_LAMBDA(int i) { data[i] = val; });
        })
    }
    else
        throw std::runtime_error("Vector::fill NO ALTERNATIVE");
}



} // namespace PHARE::core::gpu


namespace PHARE::core::mem
{

template<auto alloc_mode, typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (any_in(AllocatorMode::GPU_UNIFIED, alloc_mode))
        PHARE::core::gpu::copy(dst, src, size);
    else
        std::copy(src, src + size, dst);
}

template<auto alloc_mode0, auto alloc_mode1, typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (any_in(AllocatorMode::GPU_UNIFIED, alloc_mode0, alloc_mode1))
        PHARE::core::gpu::copy(dst, src, size);
    else
        std::copy(src, src + size, dst);
}

} // namespace PHARE::core::mem




#endif /* PHARE_CORE_UTILITIES_MEMORY_HPP */
