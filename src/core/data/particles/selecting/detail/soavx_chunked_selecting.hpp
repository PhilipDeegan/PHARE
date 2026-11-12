#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOAVX_CHUNKED_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOAVX_CHUNKED_SELECTING

// SoAVXChunked: the particles of a box are the offsets' ranges of its cells (as AoSMapped's
//  CellMap export), each particle loaded from the columns. GPU_UNIFIED reads unified memory on
//  the host

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"

namespace PHARE::core
{

template<auto alloc_mde>
struct ParticlesSelector<LayoutMode::SoAVXChunked, alloc_mde>
{
    auto constexpr static layout_mode = LayoutMode::SoAVXChunked;
    auto constexpr static alloc_mode  = alloc_mde;

    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename DstParticles, typename box_t>
    static void select(SrcParticles const& src, DstParticles& dst, box_t const& box)
    {
        src.for_each_range(box, [&](std::size_t const start, std::size_t const size) {
            for (std::size_t i = start; i < start + size; ++i)
                dst.emplace_back(src.copy(i));
        });
    }

    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename DstParticles, typename box_t, typename Shift>
    static void select(SrcParticles const& src, DstParticles& dst, box_t const& box,
                       Shift&& shifter)
    {
        src.for_each_range(box, [&](std::size_t const start, std::size_t const size) {
            for (std::size_t i = start; i < start + size; ++i)
                dst.emplace_back(shift_particle(src.copy(i), shifter));
        });
    }

    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename box_t>
    static std::size_t count(SrcParticles const& src, box_t const& box)
    {
        return src.count(box);
    }
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOAVX_CHUNKED_SELECTING */
