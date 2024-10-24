#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP

#include "core/vector.hpp"

#include "core/data/particles/selecting/particles_selecting.hpp"

namespace PHARE::core
{

template<typename ParticleArray_t, std::size_t impl = 0>
struct ParticleArraySelector
{
    using box_t = Box<int, ParticleArray_t::dimension>;
    using _impl = detail::ParticleArraySelector<ParticleArray_t, impl>;

    template<typename P2>
    void operator()(box_t const& box, P2& dst)
    {
        _impl{src}(dst, box);
    }

    template<typename P2, typename Transformer>
    /*= [](auto const& particle){ return modify(particle); }*/
    void operator()(box_t const& box, P2& dst, Transformer&& fn)
    {
        _impl{src}(dst, box, std::forward<Transformer>(fn));
    }

    ParticleArray_t const& src;
    // // scratch space for samrai boxes conversion
    // static inline thread_local std::vector<box_t> boxes;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
