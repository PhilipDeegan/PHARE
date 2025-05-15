// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER

#include "core/data/particles/exporting/detail/def_exporting.hpp"
#include "core/data/particles/exporting/particles_exporting.hpp"

namespace PHARE::core
{

template<typename Src, typename Dst, typename Box_t>
void move_particles(Src& src, Dst& dst, Box_t const& box)
{
    // copy to dst and delete from src
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID = "move_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Exporter{}.move_particles(src, dst, box);
}

template<typename Src, typename Dst, typename Box_t, typename Transformer>
void move_particles(Src& src, Dst& dst, Box_t const& box, Transformer fn0)
{
    // copy to dst and delete from src
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID = "move_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Exporter{}.move_particles(src, dst, box, fn0);
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
