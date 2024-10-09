#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"

namespace PHARE::core::detail
{

template<std::size_t dim, auto lm, auto sm, auto am, std::uint8_t impl>
struct ParticleArraySelectorImpl;


template<typename ParticleArray, std::uint8_t _impl>
class ParticleArraySelector
{
    auto static constexpr dimension    = ParticleArray::dim;
    auto static constexpr alloc_mode   = ParticleArray::alloc_mode;
    auto static constexpr layout_mode  = ParticleArray::layout_mode;
    auto static constexpr storage_mode = ParticleArray::storage_mode;

public:
    using dispatch
        = ParticleArraySelectorImpl<dimension, layout_mode, storage_mode, alloc_mode, _impl>;
    using box_t = Box<int, ParticleArray::dimension>;

    void operator()(box_t const& select)
    {
        if constexpr (ParticleArray::is_mapped)
            from.export_particles(select, to);
    }

    template<typename Transformer>
    void operator()(box_t const& select, Transformer&& fn)
    {
        if constexpr (ParticleArray::is_mapped)
            from.export_particles(select, to, std::forward<Transformer>(fn));
    }

    box_t domain_box;
    ParticleArray const& from;
    ParticleArray& to;
    std::size_t start = 0, end = from.size();
};


} // namespace PHARE::core::detail

// #include "particles_selecting_cpu.hpp"

// #if PHARE_HAVE_MKN_GPU
// #include "particles_selecting_gpu_mkn.hpp"
// #else // no impl (yet)
// namespace PHARE::core::detail
// {
// template<typename ParticleArray>
// class ParticleArraySelector<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;
// // template<typename ParticleArray>
// // class ParticleSelector<AllocatorMode::GPU, /*impl = */ 0, ParticleArray>;
// } // namespace PHARE::core::detail
// #endif

#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP*/
