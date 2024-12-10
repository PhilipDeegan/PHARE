
#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_DEF_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_DEF_CONVERTING


#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{

using enum LayoutMode;

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesConverter
{
    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<auto type, typename Dst, typename Src, typename GridLayout>
    Dst operator()(Src const& src, GridLayout const& layout);

    std::size_t const start;
    std::size_t const end;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_DEF_CONVERTING */
