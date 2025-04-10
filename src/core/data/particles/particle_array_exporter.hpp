#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER


#include "core/data/particles/exporting/detail/def_exporting.hpp"
#include "core/data/particles/exporting/particles_exporting.hpp"

namespace PHARE::core
{

// template<typename Src, typename Fn0, typename Fn1>
// void export_particles(Src const& src, Fn0 fn0, Fn1 fn1)
// {
//     using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

//     std::string_view constexpr static FN_ID = "export_particles,";
//     auto constexpr function_id              = join_string_views_v<FN_ID, Src::type_id>;
//     PHARE_LOG_SCOPE(1, function_id);

//     Exporter{}.operator()(src, fn0, fn1);

//     // PHARE_DEBUG_DO(assert(dst.size() == old_size + src.size());)
// }

template<typename Src, typename Dst, typename Box_t, typename Transformer>
void export_particles(Src const& src, Dst& dst, Box_t const& box, Transformer fn0)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID = "export_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Exporter{}.operator()(src, dst, box, fn0);
}

template<typename Src, typename Dst, typename Box_t, typename Refiner, typename Transformer>
void export_refined_particles( //
    Src const& src, Dst& dst, Box_t const& box, Refiner refiner, Transformer fn0)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID = "export_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Exporter{}.operator()(src, dst, box, refiner, fn0);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
