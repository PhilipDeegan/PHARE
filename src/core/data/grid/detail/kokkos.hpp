
#ifndef PHARE_CORE_GRID_DETAIL_KOKKOS_HPP
#define PHARE_CORE_GRID_DETAIL_KOKKOS_HPP

namespace PHARE::core::mkn_gpu
{
struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&... args)
    {
        Kokkos::parallel_for(
            "Loop1", box.size(), KOKKOS_LAMBDA(int const i) {
                auto b = box.begin();
                b      = b + i;
                fn(*b, args...);
            });
    }
};


} // namespace PHARE::core::mkn_gpu


#endif /* PHARE_CORE_GRID_DETAIL_KOKKOS_HPP */
