#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER


#include "core/data/particles/serializing/detail/def_serializing.hpp"
#include "core/data/particles/serializing/particles_serializing.hpp"

namespace PHARE::core
{


template<typename Src>
void serialize_particles(Src const& src, std::string const& file_name)
{
    using Serializer = ParticlesSerializer<Src::layout_mode, Src::alloc_mode>;

    Serializer{}.operator()(src, file_name);
}

template<typename Src, typename As = Src>
Src deserialize_particles(std::string const& file_name)
{
    using Deserializer = ParticlesDeserializer<Src::layout_mode, Src::alloc_mode>;

    // if As != Src - convert - todo

    return Deserializer{}.operator()(file_name);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER */
