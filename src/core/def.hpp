#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "core/def/phare_config.hpp"

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


namespace PHARE
{
template<typename T>
inline void throw_runtime_error([[maybe_unused]] T const& err) _PHARE_ALL_FN_
{
#if defined(__HIPCC__) || defined(__CUDACC__)
    PHARE_ASSERT(false);
#else
    throw std::runtime_error(err);
#endif
}

} // namespace PHARE



#endif // PHARE_CORE_DEF_HPP
