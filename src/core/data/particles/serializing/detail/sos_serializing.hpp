
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_SERIALIZING


#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/serializing/detail/def_serializing.hpp"
#include "core/utilities/types.hpp"
#include <optional>


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

// SoA
template<>
template<typename Src>
void ParticlesSerializer<SoA, CPU>::operator()(std::string const& file_name, Src const& src)
{
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<SoA, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    using enum LayoutMode;

    static_assert(Dst::alloc_mode == Src::alloc_mode);
    static_assert(layout_mode == Src::layout_mode || any_in(Src::layout_mode, AoS, AoSMapped));

    if constexpr (layout_mode == Src::layout_mode)
    {
        throw std::runtime_error("no impl");
    }
    else if (any_in(Src::layout_mode, AoS, AoSMapped))
    {
        //
    }
    else
        throw std::runtime_error("no impl");
}

// SoAMapped
template<>
template<typename Src>
void ParticlesSerializer<SoAVX, CPU>::operator()(std::string const& file_name, Src const& src)
{
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<SoAVX, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    using enum LayoutMode;

    static_assert(Dst::alloc_mode == Src::alloc_mode);
    static_assert(layout_mode == Src::layout_mode || any_in(Src::layout_mode, AoS, AoSMapped));

    if constexpr (layout_mode == Src::layout_mode)
    {
        throw std::runtime_error("no impl");
    }
    else if (any_in(Src::layout_mode, AoS, AoSMapped))
    {
        //
    }
    else
        throw std::runtime_error("no impl");


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_SERIALIZING */
