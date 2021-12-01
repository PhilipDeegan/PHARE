#ifndef PHARE_CORE_OPERATORS_H
#define PHARE_CORE_OPERATORS_H

#include "core/def.h"

namespace PHARE::core
{
template<typename T, bool atomic = false>
struct Operators
{
    void operator+=(T const& v)
    {
        if constexpr (atomic)
            __sync_fetch_and_add(&t, v);
        else
            t += v;
    }
    void operator+=(T const&& v) { (*this) += v; }

    void operator-=(T const& v)
    {
        if constexpr (atomic)
            __sync_fetch_and_sub(&t, v);
        else
            t += v;
    }
    void operator-=(T const&& v) { (*this) += v; }

    T& t;
};
} // namespace PHARE::core

#endif /*PHARE_CORE_OPERATORS_H*/
