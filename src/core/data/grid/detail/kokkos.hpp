#ifndef PHARE_CORE_GRID_DETAIL_KOKKOS_HPP
#define PHARE_CORE_GRID_DETAIL_KOKKOS_HPP

#include "core/def.hpp"
#include "core/utilities/box/box_spans.hpp"

namespace PHARE::core::gpu::impl
{

struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&&... args)
    {
        PHARE_LOG_LINE_SS(box.size());
        int size = box.size();
        Kokkos::parallel_for("evalOnBox", size, [&](int i) { // will fail with copy capture
            fn(point_at(box, i), args...);
        });
    }
};


} // namespace PHARE::core::gpu::impl


#endif /* PHARE_CORE_GRID_DETAIL_KOKKOS_HPP */
