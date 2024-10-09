#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <tuple>
#include <cassert>
#include <stdexcept>
#include <string_view>

#if !defined(PHARE_UNDEF_ASSERT)
#define PHARE_ASSERT(...) assert(__VA_ARGS__)
#else
#define PHARE_ASSERT(...)
#endif


#define NO_DISCARD [[nodiscard]]

// Macro String manip
#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)
#define PHARE_TOKEN_PASTE(x, y) x##y
#define PHARE_STR_CAT(x, y) PHARE_TOKEN_PASTE(x, y)


// Available Allocations
#if !defined(PHARE_WITH_GPU) or PHARE_HAVE_MKN_GPU == 1
#if defined(PHARE_HAVE_MKN_GPU) and PHARE_HAVE_MKN_GPU == 1
#define PHARE_WITH_GPU 1
#include "mkn/gpu.hpp"
#else
#define PHARE_WITH_GPU 0
#endif
#endif

#if PHARE_WITH_GPU

#define _PHARE_DEV_FN_ __device__
#define _PHARE_HST_FN_ __host__
#define _PHARE_ALL_FN_ _PHARE_HST_FN_ _PHARE_DEV_FN_

namespace PHARE
{
// !WARNING!
// PGIC++/NVC++ CANNOT HANDLE separate device/host functions!
#ifndef __HIPCC__
// #error
template<typename T>
inline void throw_runtime_error(T const&) __device__
{
    // gpu cannot throw
    PHARE_ASSERT(false);
}
#else
template<typename T>
inline void throw_runtime_error(T const& err) __host__
{
    throw std::runtime_error(err);
}
#endif
} // namespace PHARE

#else

#define _PHARE_DEV_FN_
#define _PHARE_HST_FN_
#define _PHARE_ALL_FN_

namespace PHARE
{
template<typename T>
inline void throw_runtime_error(T const& err)
{
    throw std::runtime_error(err);
}
} // namespace PHARE

#endif


#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#else
#define PHARE_DEBUG_DO(...)
#endif




#endif /* PHARE_CORE_DEF_HPP */
