// IWYU pragma: private, include "core/numerics/stable_dt/stable_dt.hpp"

#ifndef PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_MKN_GPU_HPP
#define PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_MKN_GPU_HPP

#include "core/def.hpp"
#include "core/operators.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <bit>

namespace PHARE::core::gpu::impl
{

struct StableDt
{
    // a cell per thread, atomic max of fn as the bits of a non negative double (ordered as
    //  unsigned integers), args captured by value (views)
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static double max_on_box(Fn const& fn, Box<T, S> const& box, Args const&... args)
    {
        using Op = Operators<unsigned long long, true, true>;
        gpu::Vec_t<unsigned long long> max(1, 0);
        auto* const m = max.data();
        kernel::launch(box.size(), [=] _PHARE_ALL_FN_() mutable {
            Op{*m}.max_return_old(std::bit_cast<unsigned long long>(
                static_cast<double>(fn(point_at(box, kernel::idx()), args...))));
        });
        return std::bit_cast<double>(max[0]);
    }
};

} // namespace PHARE::core::gpu::impl

#endif /* PHARE_CORE_NUMERICS_STABLE_DT_DETAIL_MKN_GPU_HPP */
