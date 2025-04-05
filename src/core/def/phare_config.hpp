#ifndef PHARE_CORE_DEF_PHARE_CONFIG_HPP
#define PHARE_CORE_DEF_PHARE_CONFIG_HPP


#include "core/def/detail/build_config.hpp"
#include "core/def/detail/umpire.hpp"
#include "core/def/detail/raja.hpp" // checks for umpire to know if gpu
#include "core/def/detail/mkn_avx.hpp"
#include "core/def/detail/mkn_gpu.hpp"
#include "core/def/kokkos_tools.hpp"
#include "core/def/detail/thrust.hpp"

#include <cstdint>

#if PHARE_HAVE_MKN_GPU or PHARE_HAVE_UMPIRE
#define PHARE_HAVE_GPU 1
#define PHARE_WITH_GPU(...) __VA_ARGS__
#else
#define PHARE_HAVE_GPU 0
#define PHARE_WITH_GPU(...)
#endif

#if PHARE_HAVE_GPU

#define _PHARE_DEV_FN_ __device__
#define _PHARE_HST_FN_ __host__
#define _PHARE_ALL_FN_ _PHARE_HST_FN_ _PHARE_DEV_FN_

#else // !PHARE_HAVE_GPU

#define _PHARE_DEV_FN_
#define _PHARE_HST_FN_
#define _PHARE_ALL_FN_

#endif // PHARE_HAVE_GPU


namespace PHARE
{
struct CompileOptions
{
    static constexpr bool WithUmpire = PHARE_HAVE_UMPIRE;
    static constexpr bool WithMknGpu = PHARE_HAVE_MKN_GPU;
    static constexpr bool WithMknAVX = PHARE_HAVE_MKN_AVX;
    static constexpr bool WithRAJA   = PHARE_HAVE_RAJA;
    static constexpr bool WithThrust = PHARE_HAVE_THRUST;
    static constexpr bool WithKokkos = PHARE_HAVE_KOKKOS;
    static constexpr bool WithGpu    = WithUmpire || WithMknGpu;
};


enum class AllocatorMode : std::uint16_t {
    CPU = 0,
    GPU_UNIFIED,
    // GPU, // unified for now
};

} // namespace PHARE

#endif /*PHARE_CORE_DEF_PHARE_CONFIG_HPP*/
