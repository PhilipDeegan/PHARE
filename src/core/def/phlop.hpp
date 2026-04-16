#ifndef PHARE_CORE_DEF_PHLOP_HPP
#define PHARE_CORE_DEF_PHLOP_HPP

#if !defined(PHARE_HAVE_PHLOP) && __has_include("phlop/timing/scope_timer.hpp")
#define PHARE_HAVE_PHLOP 1
#endif

#if !defined(PHARE_HAVE_PHLOP)
#define PHARE_HAVE_PHLOP 0
#endif //  PHARE_HAVE_PHLOP

#if PHARE_HAVE_PHLOP
#include "phlop/timing/scope_timer.hpp"
#define PHARE_WITH_PHLOP(...) __VA_ARGS__

#else
#define PHARE_WITH_PHLOP(...)

#endif // PHARE_HAVE_PHLOP


#endif /* PHARE_CORE_DEF_PHLOP_HPP */
