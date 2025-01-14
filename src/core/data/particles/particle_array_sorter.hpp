#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP

#include "core/data/particles/particle_array_def.hpp"
#include "core/vector.hpp"

#include "core/data/particles/sorting/particles_sorting.hpp"

namespace PHARE::core
{
template<typename ParticleArray, std::size_t impl = 0>
struct ParticleArraySorter
{
    using enum LayoutMode;
    using box_t     = Box<int, ParticleArray::dimension>;
    using sort_impl = detail::ParticleSorter<ParticleArray::alloc_mode, impl, ParticleArray>;

    bool constexpr static is_sortable()
    {
        return true; //! any_in(ParticleArray::layout_mode, AoSTS, SoATS);
    }

    bool constexpr static is_cell_sortable()
    {
        return !any_in(ParticleArray::layout_mode, AoSPC, SoAPC);
    }

    auto constexpr static sortable      = is_sortable();
    auto constexpr static cell_sortable = is_cell_sortable();

    void operator()()
    {
        if constexpr (sortable and cell_sortable)
            sort_impl{particles, box}();
    }

    void by_delta()
    {
        if constexpr (sortable)
            sort_impl{particles, box}.by_delta();
    }


    ParticleArray& particles;
    box_t box;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP */
