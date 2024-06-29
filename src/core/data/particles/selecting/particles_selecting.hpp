#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"

namespace PHARE::core::detail
{
template<auto alloc_mode, std::size_t impl, typename ParticleArray>
class ParticleArraySelector;
} // namespace PHARE::core::detail

#include "particles_selecting_cpu.hpp"

#if PHARE_HAVE_MKN_GPU
#include "particles_selecting_gpu_mkn.hpp"
#else // no impl (yet)
namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleArraySelector<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;
// template<typename ParticleArray>
// class ParticleSelector<AllocatorMode::GPU, /*impl = */ 0, ParticleArray>;
} // namespace PHARE::core::detail
#endif

#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP*/
