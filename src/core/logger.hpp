#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP

#include <format>
#include <cstdint>
#include <string>
#include <utility>

#if !defined(PHARE_LOG_LEVEL)
#define PHARE_LOG_LEVEL 2 // 0 == off
#endif

namespace PHARE
{
constexpr static std::uint8_t LOG_LEVEL = PHARE_LOG_LEVEL;
} // namespace PHARE

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#include <sstream>
#include <iostream>
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#define PHARE_LOG_LINE_SS(s) PHARE_LOG_LINE_STR((std::stringstream{} << s).str());
#else
#define PHARE_LOG_LINE_STR(str)
#define PHARE_LOG_LINE_SS(str)
#endif
#define PHARE_LOG_LINE PHARE_LOG_LINE_STR("")


#if PHARE_WITH_CALIPER

#include "caliper/cali.h"

#define PHARE_LOG_START(lvl, str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(lvl, str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE::scope_log __phare_scope##__line__(lvl, str)

#else // !PHARE_WITH_CALIPER

#include "core/utilities/logger/logger_defaults.hpp"


#endif // PHARE_WITH_CALIPER



namespace PHARE
{
struct scope_log
{
    scope_log(int&& i_, std::string&& str)
        : i{i_}
        , key{std::move(str)}
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_START(i, key.c_str());
        }
    }
    ~scope_log()
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_STOP(i, key.c_str());
        }
    }

    int i;
    std::string key;
};


struct ScopeTimer
{
    std::string key;
    std::size_t start_time = now();

    ~ScopeTimer()
    {
        std::cout << std::format("{:%Y-%m-%d-%H:%M:%S}", std::chrono::system_clock::now())
                  << " PHARE SCOPE TIMER: " << key << " time: " << (now() - start_time) / 1e6
                  << " ms" << std::endl;
    }

    std::uint64_t static now()
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::steady_clock::now().time_since_epoch())
            .count();
    }
};


} // namespace PHARE

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#define PHARE_FN_TIMER(key) PHARE::ScopeTimer __phare_scope##__line__(key)
#else
#define PHARE_FN_TIMER(key) // noop
#endif                      //

#endif /* PHARE_CORE_LOGGER_H */
