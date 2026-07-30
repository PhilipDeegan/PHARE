#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING


#include "core/data/particles/particle_array_def.hpp"

#include <SAMRAI/hier/BoxContainer.h>


namespace PHARE::amr
{

// dst_amr_box: destination patch's interior box, distinct from dst_boxes - lets the Ghost
// reserve estimate exclude particles whose offspring only ever land in the interior.
template<typename Src, typename Dst>
struct RefineArgs
{
    using Box_t = core::Box<int, Src::dimension>;

    Dst& dst;
    Src const& src;
    SAMRAI::hier::BoxContainer const& dst_boxes;
    Box_t const& dst_amr_box;
};

template<auto layout_mde, auto alloc_mde>
struct ParticlesRefiner
{
    static_assert(core::all_are<core::LayoutMode>(layout_mde));
    static_assert(core::all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<auto options, auto type, typename Src, typename Dst>
    void operator()(RefineArgs<Src, Dst>&, auto /*Refiner*/, auto /*Transformer*/);
};


} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING */
