
#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_DEF_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_DEF_APPENDING

#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesAppender
{
    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<auto type, typename Src, typename Dst, typename GridLayout>
    void operator()(Src const& src, Dst& dst, GridLayout const& layout);

    std::size_t const start;
    std::size_t const end;
};

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_DEF_APPENDING */