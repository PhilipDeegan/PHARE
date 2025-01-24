#ifndef PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"

namespace PHARE::core::detail
{
template<auto alloc_mode, std::size_t impl, typename ParticleArray>
class ParticleSorter;
} // namespace PHARE::core::detail

#include "particles_sorting_cpu.hpp"

#if PHARE_HAVE_MKN_GPU or PHARE_HAVE_KOKKOS
#include "particles_sorting_gpu_mkn.hpp"
#else // no impl (yet)
namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleSorter<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;
// template<typename ParticleArray>
// class ParticleSorter<alloc_mode::GPU, /*impl = */ 0, ParticleArray>;
} // namespace PHARE::core::detail
#endif

#endif /*PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP*/
