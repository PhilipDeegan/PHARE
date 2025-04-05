#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER

#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/converting/particles_converting.hpp"

namespace PHARE::core
{


template<typename Dst, typename Src, typename... Args>
auto static convert_particles(Src const& src, Args&&... args)
{
    using Converter
        = ParticlesConverter<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>;

    std::string_view constexpr static FN_ID = "convert_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    auto const& [layout, start, end] = std::forward_as_tuple(args...);
    return Converter{start, end}.template operator()<Dst>(src, layout);
}

template<typename Dst, typename Src, typename GridLayout>
auto static convert_particles(Src const& src, GridLayout const& layout)
{
    if constexpr (std::is_same_v<Dst, Src>)
        return src;
    else
        return convert_particles<Dst>(src, layout, 0ull, src.size());
}

template<typename Dst, typename Src, typename GridLayout>
auto static convert_particles_and_sort(Src const& src, GridLayout const& layout)
{
    std::string_view constexpr static FN_ID = "convert_particles_and_sort,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    auto out = convert_particles<Dst>(src, layout);
    sort_particles(out, layout.AMRBox());
    return out;
}



} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER */
