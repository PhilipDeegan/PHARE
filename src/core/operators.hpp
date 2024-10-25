#ifndef PHARE_CORE_OPERATORS_HPP
#define PHARE_CORE_OPERATORS_HPP

#include "core/def.hpp"

#ifndef PHARE_WITH_GPU
#define PHARE_WITH_GPU 0
#endif

// #include "hip/amd_detail/amd_hip_atomic.h"
#include <atomic>

#if PHARE_WITH_GPU
// #include <cuda_runtime.h>
#endif

namespace PHARE::core
{
template<typename T, bool atomic = false>
struct Operators
{
    T static constexpr ONE = 1;

    static_assert(not std::is_const_v<T>); // doesn't make sense

    auto static constexpr GPU = PHARE_WITH_GPU;

    void operator+=(T const& v) _PHARE_ALL_FN_
    {
        if constexpr (GPU and atomic)
        {
            atomicAdd(&t, v);
        }
        else if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp + v)) {}
        }
        else
            t += v;
    }
    void operator+=(T const&& v) _PHARE_ALL_FN_ { (*this) += v; }

    void operator-=(T const& v) _PHARE_ALL_FN_
    {
        if constexpr (GPU and atomic)
        {
            atomicSub(&t, v);
        }
        else if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp - v)) {}
        }
        else
            t -= v;
    }
    void operator-=(T const&& v) _PHARE_ALL_FN_ { (*this) += v; }

    auto increment_return_old() _PHARE_ALL_FN_ // postfix increment
    {
        if constexpr (GPU and atomic)
        {
            auto o = atomicAdd(&t, ONE);
            PHARE_ASSERT(o < t);
            return o;
        }
        else if constexpr (atomic)
        {
            throw std::runtime_error("finish");
            // auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            // T tmp          = atomic_t.load();
            // while (!atomic_t.compare_exchange_weak(tmp, tmp + v)) {}
        }
        else
        {
            T tmp = t;
            ++t;
            return tmp;
        }
    }

    T& t;
};
} // namespace PHARE::core

#endif /* PHARE_CORE_OPERATORS_HPP */
