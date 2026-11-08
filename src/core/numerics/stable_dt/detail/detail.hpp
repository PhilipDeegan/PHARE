// IWYU pragma: private, include "core/numerics/stable_dt/stable_dt.hpp"

#ifndef PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_DETAIL_HPP
#define PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_DETAIL_HPP

#include "core/def.hpp"

#if PHARE_HAVE_MKN_GPU
#include "core/numerics/stable_dt/detail/mkn_gpu.hpp"
#endif

#if PHARE_HAVE_KOKKOS
#include "core/numerics/stable_dt/detail/kokkos.hpp"
#endif

#include "core/utilities/box/box.hpp"

#include <cstddef>
#include <stdexcept>

namespace PHARE::core::gpu
{
struct StableDt
{
    // largest fn(cell, args...) over box, fn must return a non negative value
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static double max_on_box(Fn const& fn, Box<T, S> const& box, Args const&... args)
    {
        PHARE_WITH_KOKKOS({ //
            return impl::StableDt::max_on_box(fn, box, args...);
        })

        PHARE_WITH_MKN_GPU({ //
            return impl::StableDt::max_on_box(fn, box, args...);
        })

        throw std::runtime_error("gpu::StableDt::max_on_box: no gpu backend");
    }
};

} // namespace PHARE::core::gpu


#endif /* PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_DETAIL_HPP */
