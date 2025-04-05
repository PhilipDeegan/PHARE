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

#if PHARE_HAVE_GPU

#include "particles_partitioning_gpu.hpp"

#else // no impl (yet)

namespace PHARE::core::detail
{

template<typename ParticleArray> // force compile error showing no impl
class ParticleArrayPartitioner<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;

} // namespace PHARE::core::detail

#endif

#endif /*PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP*/
