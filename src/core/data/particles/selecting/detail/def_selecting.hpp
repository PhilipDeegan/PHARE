
#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_DEF_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_DEF_SELECTING

#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesSelector
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename SrcParticles, typename DstParticles, typename box_t>
    static void select(SrcParticles const&, DstParticles&, box_t const&);

    template<typename SrcParticles,typename DstParticles, typename box_t, typename Shift>
    static void select(SrcParticles const&, DstParticles&, box_t const&, Shift&&);
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_DEF_SELECTING */
