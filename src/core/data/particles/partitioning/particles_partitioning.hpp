#ifndef PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"

#include "core/utilities/partitionner/partitionner.hpp"

namespace PHARE::core::detail
{
template<auto alloc_mode, std::size_t impl, typename ParticleArray>
class ParticleArrayPartitioner;
} // namespace PHARE::core::detail

#include "particles_partitioning_cpu.hpp"

#if PHARE_HAVE_MKN_GPU
#include "particles_partitioning_gpu_mkn.hpp"
#else // no impl (yet)
namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleArrayPartitioner<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;
// template<typename ParticleArray>
// class ParticlePartitioner<AllocatorMode::GPU, /*impl = */ 0, ParticleArray>;
} // namespace PHARE::core::detail
#endif

#endif /*PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP*/
