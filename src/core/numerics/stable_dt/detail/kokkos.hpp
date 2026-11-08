// IWYU pragma: private, include "core/numerics/stable_dt/stable_dt.hpp"

#ifndef PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_KOKKOS_HPP
#define PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_KOKKOS_HPP

#include "core/def.hpp"
#include "core/utilities/box/box_spans.hpp"

namespace PHARE::core::gpu::impl
{

struct StableDt
{
    // args captured by value (views)
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static double max_on_box(Fn const& fn, Box<T, S> const& box, Args const&... args)
    {
        double max = 0;
        Kokkos::parallel_reduce(
            "max_on_box", box.size(),
            KOKKOS_LAMBDA(std::size_t const i, double& m) {
                double const v = fn(point_at(box, i), args...);
                if (v > m)
                    m = v;
            },
            Kokkos::Max<double>(max));
        return max;
    }
};

} // namespace PHARE::core::gpu::impl

#endif /* PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_KOKKOS_HPP */
