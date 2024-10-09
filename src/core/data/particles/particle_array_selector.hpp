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

    void operator()(box_t const& box) { _impl{particles, box}(); }

    template<typename P2>
    void operator()(box_t const& box, P2& into)
    {
        _impl{particles}(box, into);
    }

    template<typename P2, typename Transformer>
    // Transformer == [](auto const& particle){ return modify(particle); }
    void operator()(box_t const& box, P2& into, Transformer&& fn)
    {
        _impl{particles}(box, into, std::forward<Transformer>(fn));
    }

    ParticleArray_t const& particles;

    // scratch space for samrai boxes conversion
    static inline thread_local std::vector<box_t> boxes;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
