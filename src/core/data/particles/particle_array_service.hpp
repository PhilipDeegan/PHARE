#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP

#include <array>
#include <cstdint>

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"
#include "core/vector.hpp"

namespace PHARE::core
{


struct ParticleArrayService
{
    // optionally implemented on arrays
    template<std::uint8_t Phase, auto type, typename ParticleArray_t>
    auto _PHARE_ALL_FN_ static sync(ParticleArray_t& p) -> decltype(p.sync(), void())
    {
        PHARE_LOG_SCOPE(2, "ParticleArrayService::sync");

        p.template sync<Phase, type>();
    }

    template<auto type, typename ParticleArray_t>
    auto static reserve_ppc_in(ParticleArray_t& p, std::size_t ppc)
        -> decltype(p.template reserve_ppc<type>(0ull), void())
    {
        PHARE_LOG_SCOPE(2, "ParticleArrayService::reserve_ppc_in");

        using enum LayoutMode;
        if constexpr (any_in(ParticleArray_t::layout_mode, AoSPC, SoAPC, AoSTS, SoATS))
            p.template reserve_ppc<type>(ppc);
    }




    // noops below


    template<std::uint8_t Phase, auto type, typename... Args>
    void _PHARE_ALL_FN_ static sync(Args&&...)
    {
        PHARE_LOG_SCOPE(2, "ParticleArrayService::sync noop");
    }


    template<auto type, typename... Args>
    void static reserve_ppc_in(Args&&...)
    {
        PHARE_LOG_SCOPE(2, "ParticleArrayService::reserve_ppc_in noop");
    }


    //
};



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP*/
