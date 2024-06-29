

#include "python3/cpp_simulator.hpp"
#include "core/data/particles/particle_array_detail.hpp"

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_sim_2_1_4, m)
{
    using dim                     = std::integral_constant<std::size_t, 2>;
    using interp                  = std::integral_constant<std::size_t, 1>;
    using nbRefinePart            = std::integral_constant<std::size_t, 4>;
    auto constexpr layout         = core::LayoutMode::AoSMapped;
    auto constexpr allocator_mode = AllocatorMode::CPU;

    declare_essential(m);
    declare_sim<dim, interp, nbRefinePart, layout, allocator_mode>(m);
}
} // namespace PHARE::pydata
