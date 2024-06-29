#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_CPU_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_CPU_HPP

// no includes, is included

namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleArraySelector<AllocatorMode::CPU, /*impl = */ 0, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    void operator()(box_t const& select)
    {
        if constexpr (ParticleArray::is_mapped)
            from.export_particles(select, to);
    }

    template<typename Transformation>
    void operator()(box_t const& select, Transformation&& transformer)
    {
        if constexpr (ParticleArray::is_mapped)
            from.export_particles(select, to, std::forward<Transformation>(transformer));
    }

    box_t domain_box;
    ParticleArray const& from;
    ParticleArray& to;
    std::size_t start = 0, end = from.size();
};



} // namespace PHARE::core::detail

#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_CPU_HPP*/
