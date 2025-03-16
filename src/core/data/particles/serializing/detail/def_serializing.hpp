
#ifndef PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING

#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesSerializer
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Src>
    void operator()(std::string const& file_name, Src const& src);
};



template<auto layout_mde, auto alloc_mde>
struct ParticlesDeserializer
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Dst, typename Src = Dst>
    void operator()(std::string const& file_name, Dst& dst);
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING */
