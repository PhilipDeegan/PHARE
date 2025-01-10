#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP

#include "core/data/particles/selecting/particles_selecting.hpp"

namespace PHARE::core
{

template<typename Src, typename Dst, typename box_t>
void select_particles(Src const& src, Dst& dst, box_t const& box)
{
    static_assert(std::is_same_v<Src, Dst>);
    using Selector = ParticlesSelector<Src::layout_mode, Src::alloc_mode>;

    PHARE_LOG_LINE_SS("");
    std::string_view constexpr static FN_ID = "select_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Selector::select(src, dst, box);
}


template<typename Src, typename Dst, typename box_t, typename Transformer>
void select_particles(Src const& src, Dst& dst, box_t const& box, Transformer&& transformer)
{
    static_assert(std::is_same_v<Src, Dst>);
    using Selector = ParticlesSelector<Src::layout_mode, Src::alloc_mode>;

    PHARE_LOG_LINE_SS("");
    std::string_view constexpr static FN_ID = "select_transformed_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Selector::select(src, dst, box, transformer);
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
