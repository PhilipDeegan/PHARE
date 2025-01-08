
#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING

#include "core/utilities/memory.hpp"
#include "core/data/particles/converting/detail/def_converting.hpp"

namespace PHARE::core
{

using enum LayoutMode;

template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, AllocatorMode::CPU, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
        std::copy(tile().begin(), tile().end(), std::back_inserter(out));

    return out;
}


template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, AllocatorMode::GPU_UNIFIED, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.resize(src.size());
    std::size_t offset = 0;
    for (auto const& tile : src())
    {
        auto& particles = src(src.local_cell(tile.lower));
        gpu::copy(out.data() + offset, particles.data(), particles.size());
        offset += particles.size();
    }
    return out;
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING */
