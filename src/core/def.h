#ifndef PHARE_CORE_DEF_H
#define PHARE_CORE_DEF_H

#include <cassert>
#include <stdexcept>
#include <string_view>

#if !defined(PHARE_WITH_GPU)
#if defined(HAVE_RAJA) and defined(HAVE_UMPIRE)
#define PHARE_WITH_GPU 1
#endif
#endif

#if defined(PHARE_WITH_GPU)

#define _PHARE_GPU_FN_DEV_ __device__
#define _PHARE_GPU_FN_HST_ __host__
#define _PHARE_FN_SIG_ _PHARE_GPU_FN_HST_ _PHARE_GPU_FN_DEV_

namespace PHARE
{
/*
template<typename T>
inline void throw_runtime_error(char const*const) __device__
{
    // gpu cannot throw
    assert(false);
}*/

inline void throw_runtime_error(char const*const err) __host__
{
    throw std::runtime_error(err);
}
} // namespace PHARE

#else

#define _PHARE_GPU_FN_DEV_
#define _PHARE_GPU_FN_HST_
#define _PHARE_FN_SIG_

namespace PHARE
{
template<typename T>
inline void throw_runtime_error(T const& err)
{
    throw std::runtime_error(err);
}
} // namespace PHARE


#endif


#endif /*PHARE_CORE_DEF_H*/
