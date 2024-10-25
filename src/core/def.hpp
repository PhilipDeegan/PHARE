#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <cassert>
#include <stdexcept>
#include <type_traits>


#define NO_DISCARD [[nodiscard]]

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#else
#define PHARE_DEBUG_DO(...)
#endif


#if !defined(PHARE_UNDEF_ASSERT)
//  Cuda can fail to compile with assertions
//  I've seen a github issue, will ref
#define PHARE_ASSERT(...) assert(__VA_ARGS__)
#else
#define PHARE_ASSERT(...)
#endif


#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)
#define PHARE_TOKEN_PASTE(x, y) x##y
#define PHARE_STR_CAT(x, y) PHARE_TOKEN_PASTE(x, y)


namespace PHARE::core
{

NO_DISCARD bool isUsable(auto const&... args)
{
    auto check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg != nullptr;
        else
            return arg.isUsable();
    };
    return (check(args) && ...);
}


NO_DISCARD bool isSettable(auto const&... args)
{
    auto check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg == nullptr;
        else
            return arg.isSettable();
    };
    return (check(args) && ...);
}

} // namespace PHARE::core



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
// #if !defined(__HIPCC__) || !defined(__CUDACC__)
// #error
template<typename T>
inline void throw_runtime_error(T const&) __device__
{
    // gpu cannot throw
    PHARE_ASSERT(false);
}

// #else // PHARE_WITH_GPU

template<typename T>
inline void throw_runtime_error(T const& err) __host__
{
    throw std::runtime_error(err);
}

// #endif

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

#endif // PHARE_WITH_GPU


#endif // PHARE_CORE_DEF_HPP
