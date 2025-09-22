
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_DEF_EXPORTING
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_DEF_EXPORTING

#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesExporter
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Src, typename Dst, typename Box_t>
    void move_particles(Src& src, Dst& dst, Box_t const& box, std::size_t const growby = 0);


    template<bool in = true, typename Src, std::size_t dim>
    void delete_particles(Src& src, Box<int, dim> const& box);

    template<bool in = true, typename Src, typename Boxes>
    void delete_particles(Src& src, Boxes const& boxes);

    // template<typename Src, typename Dst, typename Box_t, typename Refiner, typename Transformer>
    // void operator()(Src const&, Dst&, Box_t const&, Refiner, Transformer);
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_DEF_EXPORTING */
