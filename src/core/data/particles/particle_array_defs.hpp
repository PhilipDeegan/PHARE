#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEFS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEFS_HPP


#if !defined(PHARE_USE_CELLMAP)
#define PHARE_USE_CELLMAP 0
#endif /*PHARE_USE_CELLMAP*/

namespace PHARE::core
{
struct ParticleArrayDetails
{
    bool static constexpr use_cellmap = PHARE_USE_CELLMAP;
};
} // namespace PHARE::core

#endif
