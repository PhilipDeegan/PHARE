#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP


#include "core/data/particles/selecting/particles_selecting.hpp"

namespace PHARE::core
{

template<typename ParticleArray_t, std::size_t impl = 0>
struct ParticleArraySelector
{
    auto static constexpr dimension    = ParticleArray_t::dimension;
    auto static constexpr alloc_mode   = ParticleArray_t::alloc_mode;
    auto static constexpr layout_mode  = ParticleArray_t::layout_mode;
    auto static constexpr storage_mode = ParticleArray_t::storage_mode;

    using box_t = Box<int, dimension>;
    using _impl = detail::ParticleArraySelector<alloc_mode, layout_mode, impl, ParticleArray_t>;

    template<typename P2>
    void operator()(box_t const& box, P2& dst)
    {
        _impl::select(src, dst, box);
    }

    template<typename P2, typename Transformer>
    /* Transformer = [](auto const& particle){ return modify(particle); }*/
    void operator()(box_t const& box, P2& dst, Transformer&& fn)
    {
        _impl::select(src, dst, box, std::forward<Transformer>(fn));
    }

    ParticleArray_t const& src;
    std::size_t start = 0, end = src.size();

    // // scratch space for samrai boxes conversion
    // static inline thread_local std::vector<box_t> boxes;
};


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
