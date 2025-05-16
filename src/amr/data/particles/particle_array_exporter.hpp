// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER

#include "amr/data/particles/exporting/detail/def_exporting.hpp"
#include "amr/data/particles/exporting/particles_exporting.hpp"

namespace PHARE::amr
{


template<auto type, typename Src, typename Dst, typename Box_t, typename Refiner,
         typename Transformer>
void export_refined_particles( //
    Src const& src, Dst& dst, Box_t const& box, Refiner refiner, Transformer fn0)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID = "export_particles,";
    auto constexpr function_id              = core::join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    Exporter{}.template operator()<type>(src, dst, box, refiner, fn0);
}

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
