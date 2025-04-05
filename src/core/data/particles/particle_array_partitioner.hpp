#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP

#include "core/vector.hpp"

#include "core/data/particles/partitioning/particles_partitioning.hpp"

namespace PHARE::core
{
template<typename ParticleArray, std::size_t impl = 0>
struct ParticleArrayPartitioner
{
    using box_t = Box<int, ParticleArray::dimension>;
    using _impl = detail::ParticleArrayPartitioner<ParticleArray::alloc_mode, impl, ParticleArray>;


    template<typename... Args>
    auto operator()(Args&&... args)
    {
        return _impl{particles}(std::forward<Args>(args)...);
    }

    template<typename... Args>
    auto notIn(Args&&... args)
    {
        return _impl{particles}.notIn(std::forward<Args>(args)...);
    }

    ParticleArray& particles;
    // box_t box;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP */
