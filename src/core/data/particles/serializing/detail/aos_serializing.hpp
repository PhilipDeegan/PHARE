
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_SERIALIZING


#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/serializing/detail/def_serializing.hpp"
#include <optional>


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

// AoS
template<>
template<typename Src>
void ParticlesSerializer<AoS, CPU>::operator()(std::string const& file_name, Src const& src)
{
    using Particle_t = typename Src::value_type;
    std::ofstream f{file_name, std::ios::binary};
    f.write(reinterpret_cast<char const*>(src.vector().data()), src.size() * sizeof(Particle_t));
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoS, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    static_assert(Dst::alloc_mode == Src::alloc_mode);
    static_assert(layout_mode == Src::layout_mode || any_in(Src::layout_mode, AoS, AoSMapped));

    using Particle_t = typename Dst::value_type;

    std::ifstream f{file_name, std::ios::binary};

    // Stop eating new lines in binary mode
    f.unsetf(std::ios::skipws);

    // get its size:
    std::streampos fileSize;
    f.seekg(0, std::ios::end);
    fileSize = f.tellg();
    f.seekg(0, std::ios::beg);

    auto const curr_size = dst.size();
    dst.resize(curr_size + (fileSize / sizeof(Particle_t)));

    // read the data:
    f.read(reinterpret_cast<char*>(dst.vector().data() + curr_size),
           dst.size() * sizeof(Particle_t));
}

// AoSMapped
template<>
template<typename Src>
void ParticlesSerializer<AoSMapped, CPU>::operator()(std::string const& file_name, Src const& src)
{
    ParticlesSerializer<AoS, CPU>{}(file_name, src);
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoSMapped, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    ParticlesDeserializer<AoS, CPU>{}(file_name, dst);
    dst.remap();
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_SERIALIZING */
