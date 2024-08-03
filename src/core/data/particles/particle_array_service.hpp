#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP

#include <array>
#include <cstdint>

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{


struct ParticleArrayService
{
    template<typename P_t0, typename P_t1>
    void static copy_ghost_into_domain(P_t0 const& ghost, P_t1& domain)
    {
        static_assert(P_t0::layout_mode == LayoutMode::AoSPC);
        static_assert(P_t0::layout_mode == P_t1::layout_mode);

        domain.insert_domain_from(ghost);
    }

    template<typename P_t0, typename P_t1>
    void static copy_into(P_t0 const& src, P_t1& dst)
    {
        dst.insert(src);
    }

    // optionally implemented on arrays
    template<std::uint8_t Phase, auto type, typename ParticleArray_t>
    auto _PHARE_ALL_FN_ static sync(ParticleArray_t& p) -> decltype(p.sync(), void())
    {
        p.template sync<Phase, type>();
    }

    template<auto type, typename ParticleArray_t>
    auto static reserve_ppc_in(ParticleArray_t& p,
                               std::size_t ppc) -> decltype(p.template reserve_ppc<type>(0ull),
                                                            void())
    {
        if constexpr (ParticleArray_t::layout_mode == LayoutMode::AoSPC)
            p.template reserve_ppc<type>(ppc);
    }



    // noops below
    template<std::uint8_t Phase, auto type, typename... Args>
    void _PHARE_ALL_FN_ static sync(Args&&...)
    {
    }


    template<auto type, typename... Args>
    void static reserve_ppc_in(Args&&...)
    {
    }
};



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP*/
