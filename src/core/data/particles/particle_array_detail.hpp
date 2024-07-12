#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP

#include <cstdint>

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_aos.hpp"
#include "core/data/particles/arrays/particle_array_aos_pc.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soa_pc.hpp"

namespace PHARE::core
{

template<auto layout_mode, std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver;


template<std::size_t dim,                            //
         auto layout_mode_  = LayoutMode::AoSMapped, //
         auto storage_mode_ = StorageMode::VECTOR,   //
         auto alloc_mode_   = AllocatorMode::CPU,    //
         std::uint8_t impl_ = 0>
class ParticleArrayInternals
{
    static_assert(std::is_same_v<decltype(alloc_mode_), AllocatorMode>);
    static_assert(std::is_same_v<decltype(layout_mode_), LayoutMode>);
    static_assert(std::is_same_v<decltype(storage_mode_), StorageMode>);

    auto static constexpr resolve_t()
    {
        return ParticleArrayLayoutResolver<layout_mode, dim, storage_mode, alloc_mode,
                                           impl>::resolve_t();
    }

    struct strings
    {
        std::string_view static constexpr alloc_mode   = magic_enum::enum_name(alloc_mode_);
        std::string_view static constexpr layout_mode  = magic_enum::enum_name(layout_mode_);
        std::string_view static constexpr storage_mode = magic_enum::enum_name(storage_mode_);
        std::string_view static constexpr cma          = ",";
        std::string_view static constexpr _dim         = to_string_view_v<std::size_t, dim>;
        std::string_view static constexpr _impl        = to_string_view_v<std::uint8_t, impl_>;

        std::string_view static constexpr type_id
            = join_string_views_v<_dim, cma, layout_mode, cma, alloc_mode, cma, storage_mode,cma, _impl>;
    };

public:
    auto static constexpr impl         = impl_;
    auto static constexpr layout_mode  = layout_mode_;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = storage_mode_;
    auto static constexpr is_mapped // update accordingly
        = layout_mode == LayoutMode::AoSMapped /*|| layout_mode == LayoutMode::SoAMapped*/;

    auto static constexpr type_id = strings::type_id;

    using value_type = std::decay_t<decltype(*resolve_t())>;
};


namespace
{
    template<typename T>
    auto static constexpr _as_nullptr_()
    {
        return T{nullptr};
    }
} // namespace

template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::AoS, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<AoSParticles<AoSVector<dim, alloc_mode>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<AoSParticles<AoSSpan<dim>>*>();
    }
};




template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::AoSMapped, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<AoSMappedParticles<AoSParticles<AoSVector<dim, alloc_mode>>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<AoSMappedParticles<AoSParticles<AoSSpan<dim>>>*>();
    }
};




template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::SoA, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<SoAParticles<SoAVector<dim, alloc_mode>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<SoAParticles<SoASpan<dim /*, is_const*/>>*>();
    }
};


template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::AoSPC, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<AoSPCParticles<AoSPCVector<dim, alloc_mode, impl>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<AoSPCParticles<AoSPCSpan<dim, alloc_mode, impl>>*>();
    }
};




} // namespace PHARE::core



#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP*/
