#ifndef PHARE_CORE_GRID_DETAIL_DETAIL_HPP
#define PHARE_CORE_GRID_DETAIL_DETAIL_HPP

// one kernel backend: mkn::gpu launches when present (also on umpire memory), else raja, else kokkos
#if PHARE_HAVE_MKN_GPU
#include "core/data/grid/detail/mkn_gpu.hpp"
#elif PHARE_HAVE_RAJA
#include "core/data/grid/detail/raja.hpp"
#elif PHARE_HAVE_KOKKOS
#include "core/data/grid/detail/kokkos.hpp"
#endif

#include "core/utilities/box/box.hpp"

#include <cstddef>

namespace PHARE::core::gpu
{
struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&... args)
    {
        PHARE_LOG_LINE_SS(box.size());


#if PHARE_HAVE_MKN_GPU || PHARE_HAVE_RAJA || PHARE_HAVE_KOKKOS
        impl::GridLayout::evalOnBox(fn, box, args...);
#endif
    }
};


} // namespace PHARE::core::gpu



#endif /* PHARE_CORE_GRID_DETAIL_DETAIL_HPP */
