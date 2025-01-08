#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP

#include <cstdint>

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_aos.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"
#include "core/data/particles/arrays/particle_array_pc.hpp"
#include "core/data/particles/arrays/particle_array_ts.hpp"



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

        auto static constexpr type_id = join_string_views_v<_dim, cma, layout_mode, cma, alloc_mode,
                                                            cma, storage_mode, cma, _impl>;
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

template<std::size_t dim, typename internals = ParticleArrayInternals<dim>>
class ParticleArray;

template<std::size_t d,     // dim
         auto l,            // layout
         auto s,            // storage
         auto a,            // alloc
         std::uint8_t i = 0 // impl
         >
using ParticleArrayInternals_vt = typename ParticleArrayInternals<d, l, s, a, i>::value_type;



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
            return _as_nullptr_<AoSParticles<AoSSpan<dim, alloc_mode>>*>();
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
            return _as_nullptr_<AoSMappedParticles<AoSParticles<AoSSpan<dim, alloc_mode>>>*>();
    }
};




template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::AoSPC, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        using Inner
            = ParticleArrayInternals_vt<dim, LayoutMode::AoS, storage_mode, alloc_mode, impl>;

        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<PerCellParticles<PerCellVector<Inner, impl>>*>();

        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<PerCellParticles<PerCellSpan<Inner, impl>>*>();
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
            return _as_nullptr_<SoAParticles<SoASpan<dim, alloc_mode>>*>();
    }
};



template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::SoAPC, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        using Inner
            = ParticleArrayInternals_vt<dim, LayoutMode::SoA, storage_mode, alloc_mode, impl>;

        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<PerCellParticles<PerCellVector<Inner, impl>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<PerCellParticles<PerCellSpan<Inner, impl>>*>();
    }
};


template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::AoSTS, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<
            dim, ParticleArrayInternals<dim, LayoutMode::AoS, storage_mode, alloc_mode, impl>>;

        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, impl>>*>();

        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, impl>>*>();
    }
};


template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::SoATS, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<
            dim, ParticleArrayInternals<dim, LayoutMode::SoA, storage_mode, alloc_mode, impl>>;

        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, impl>>*>();

        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, impl>>*>();
    }
};


template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::SoAVX, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<SoAVXParticles<SoAVXVector<dim, alloc_mode>>*>();
        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<SoAVXParticles<SoAVXSpan<dim, alloc_mode>>*>();
    }
};
template<std::size_t dim, auto storage_mode, auto alloc_mode, std::uint8_t impl>
struct ParticleArrayLayoutResolver<LayoutMode::SoAVXTS, dim, storage_mode, alloc_mode, impl>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<
            dim, ParticleArrayInternals<dim, LayoutMode::SoAVX, storage_mode, alloc_mode, impl>>;

        if constexpr (storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, impl>>*>();

        if constexpr (storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, impl>>*>();
    }
};


} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP*/
