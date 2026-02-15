#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP

#include "core/data/particles/particle_array_def.hpp"

// Mutators
#include "core/data/particles/arrays/particle_array_pc.hpp"
#include "core/data/particles/arrays/particle_array_ts.hpp"

// Impls
#include "core/data/particles/arrays/particle_array_aos.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"
#include "core/def/phare_config.hpp"

#include <cstdint>


namespace PHARE::core
{


struct ParticleArrayOptions
{
    std::size_t dim          = 1;
    LayoutMode layout_mode   = LayoutMode::AoSMapped;
    StorageMode storage_mode = StorageMode::VECTOR;
    AllocatorMode alloc_mode = AllocatorMode::CPU;
    bool _const_             = 0; // sometimes needed
    std::uint8_t impl        = 2; // TORM;

    auto constexpr with_layout(LayoutMode const lm) const
    {
        auto copy        = *this;
        copy.layout_mode = lm;
        return copy;
    }
    auto constexpr with_storage(StorageMode const sm) const
    {
        auto copy         = *this;
        copy.storage_mode = sm;
        return copy;
    }
    auto constexpr with_alloc(AllocatorMode const am) const
    {
        auto copy       = *this;
        copy.alloc_mode = am;
        return copy;
    }
};




template<auto layout_mode, auto opts>
struct ParticleArrayLayoutResolver;


template<auto opts>
class ParticleArrayInternals
{
    static_assert(std::is_same_v<decltype(opts.alloc_mode), AllocatorMode>);
    static_assert(std::is_same_v<decltype(opts.layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(opts.storage_mode), StorageMode>);

    auto static constexpr resolve_t()
    {
        return ParticleArrayLayoutResolver<opts.layout_mode, opts>::resolve_t();
    }

    struct strings
    {
        std::string_view static constexpr alloc_mode   = magic_enum::enum_name(opts.alloc_mode);
        std::string_view static constexpr layout_mode  = magic_enum::enum_name(opts.layout_mode);
        std::string_view static constexpr storage_mode = magic_enum::enum_name(opts.storage_mode);
        std::string_view static constexpr cma          = ",";
        std::string_view static constexpr _dim         = to_string_view_v<std::size_t, opts.dim>;

        auto static constexpr type_id
            = join_string_views_v<_dim, cma, layout_mode, cma, alloc_mode, cma, storage_mode, cma>;
    };

public:
    auto static constexpr impl = opts.impl;
    auto static constexpr is_mapped // update accordingly
        = opts.layout_mode
          == LayoutMode::AoSMapped /*|| layout_mode == LayoutMode::SoAMapped*/; // torm?

    auto static constexpr type_id = strings::type_id;

    using value_type = std::decay_t<decltype(*resolve_t())>;
};




template<auto opts, typename internals = ParticleArrayInternals<opts>>
class ParticleArray;

template<auto opts>
using ParticleArrayInternals_vt = ParticleArrayInternals<opts>::value_type;



namespace
{
    template<typename T>
    auto static constexpr _as_nullptr_()
    {
        return T{nullptr};
    }
} // namespace



template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::AoS, opts>
{
    auto static constexpr resolve_t()
    {
        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<AoSParticles<AoSVector<opts.dim, opts.alloc_mode>>*>();
        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<AoSParticles<AoSSpan<opts.dim, opts.alloc_mode>>*>();
    }
};




template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::AoSMapped, opts>
{
    auto static constexpr resolve_t()
    {
        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<
                AoSMappedParticles<AoSParticles<AoSMappedVector<opts.dim, opts.alloc_mode>>>*>();
        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<
                AoSMappedParticles<AoSParticles<AoSMappedSpan<opts.dim, opts.alloc_mode>>>*>();
    }
};




template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::AoSPC, opts>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<opts.with_layout(LayoutMode::AoS)>;

        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<PerCellParticles<PerCellVector<Inner, opts.impl>>*>();

        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<PerCellParticles<PerCellSpan<Inner, opts.impl>>*>();
    }
};




template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::SoA, opts>
{
    auto static constexpr resolve_t()
    {
        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<SoAParticles<SoAVector<opts.dim, opts.alloc_mode>>*>();
        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<SoAParticles<SoASpan<opts.dim, opts.alloc_mode, opts._const_>>*>();
    }
};



template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::SoAPC, opts>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<opts.with_layout(LayoutMode::SoA)>;

        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<PerCellParticles<PerCellVector<Inner, opts.impl>>*>();
        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<PerCellParticles<PerCellSpan<Inner, opts.impl>>*>();
    }
};


template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::AoSTS, opts>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<opts.with_layout(LayoutMode::AoS)>;

        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, opts.impl>>*>();

        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, opts.impl>>*>();
    }
};



template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::SoATS, opts>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<opts.with_layout(LayoutMode::SoA)>;

        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, opts.impl>>*>();

        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, opts.impl>>*>();
    }
};


template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::SoAVX, opts>
{
    auto static constexpr resolve_t()
    {
        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<SoAVXParticles<SoAVXVector<opts.dim, opts.alloc_mode>>*>();
        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<SoAVXParticles<SoAVXSpan<opts.dim, opts.alloc_mode>>*>();
    }
};
template<auto opts>
struct ParticleArrayLayoutResolver<LayoutMode::SoAVXTS, opts>
{
    auto static constexpr resolve_t()
    {
        using Inner = ParticleArray<opts.with_layout(LayoutMode::SoAVX)>;

        if constexpr (opts.storage_mode == StorageMode::VECTOR)
            return _as_nullptr_<TileSetParticles<TileSetVector<Inner, opts.impl>>*>();

        if constexpr (opts.storage_mode == StorageMode::SPAN)
            return _as_nullptr_<TileSetParticles<TileSetSpan<Inner, opts.impl>>*>();
    }
};



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP*/
