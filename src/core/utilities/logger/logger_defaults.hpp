#ifndef PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP
#define PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP

#include "core/def/phlop.hpp"
#include "core/def/kokkos_tools.hpp"


#if PHARE_HAVE_PHLOP && PHARE_HAVE_KOKKOS_TOOLS
#define PHARE_SCOPE_TIMER(str)                                                                     \
    PHLOP_SCOPE_TIMER(str);                                                                        \
    PHARE_STR_CAT(Kokkos::Profiling::ScopedRegion __phare_scope_, __line__)(str)
#endif // PHARE_HAVE_PHLOP && PHARE_HAVE_KOKKOS_TOOLS


#if !defined(PHARE_SCOPE_TIMER) && PHARE_HAVE_PHLOP
#define PHARE_SCOPE_TIMER PHLOP_SCOPE_TIMER
#endif // PHARE_HAVE_PHLOP


#if !defined(PHARE_SCOPE_TIMER) && PHARE_HAVE_KOKKOS_TOOLS
#define PHARE_SCOPE_TIMER PHARE_STR_CAT(Kokkos::Profiling::ScopedRegion __phare_scope_, __line__)
#endif // PHARE_HAVE_PHLOP


#ifndef PHARE_SCOPE_TIMER
#define PHARE_SCOPE_TIMER(str)
#endif // PHARE_SCOPE_TIMER


#if PHARE_LOG_LEVEL >= 1
#define PHARE_LOG_SCOPE_1(str) PHARE_SCOPE_TIMER(str)
#else
#define PHARE_LOG_SCOPE_1(str)
#endif // LOG_LEVEL >= 1


#if PHARE_LOG_LEVEL >= 2
#define PHARE_LOG_SCOPE_2(str) PHARE_SCOPE_TIMER(str)
#else
#define PHARE_LOG_SCOPE_2(str)
#endif // LOG_LEVEL >= 2


#if PHARE_LOG_LEVEL == 3
#define PHARE_LOG_SCOPE_3(str) PHARE_SCOPE_TIMER(str)
#else
#define PHARE_LOG_SCOPE_3(str)
#endif // LOG_LEVEL == 3

#define PHARE_LOG_START(lvl, str)
#define PHARE_LOG_STOP(lvl, str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE_STR_CAT(PHARE_LOG_SCOPE_, lvl)(str)


#endif /* PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP */
