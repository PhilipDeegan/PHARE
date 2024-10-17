#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_GPU_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_GPU_HPP

#include "core/utilities/box/box.hpp"

#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>


namespace PHARE::core::detail
{

template<typename ParticleArray>
class ParticleArraySelector<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    void operator()(std::uint64_t l, std::uint64_t r)
    {
        // PHARE_LOG_LINE_STR("");
        if constexpr (ParticleArray::layout_mode == LayoutMode::AoS)
        {
            // auto ps = particles.view();
            // thrust::sort(thrust::device, ps.begin() + l, ps.begin() + r,
            //              [cf = cell_flattener] _PHARE_ALL_FN_(auto const& a, auto const& b) {
            //                  return cf(a.iCell()) < cf(b.iCell());
            //              });
        }
        else
        {
            // auto it = SoAIteratorAdaptor::make(particles);
            // thrust::sort(thrust::device, it + l, it + r,
            //              [cf = cell_flattener] _PHARE_ALL_FN_(auto const& a, auto const& b) {
            //                  return cf(SoAIteratorAdaptor::iCell(a))
            //                         < cf(SoAIteratorAdaptor::iCell(b));
            //              });
        }
    }

    void operator()() { (*this)(0, particles.size()); }



    ParticleArray& particles;
    box_t& domain_box;
    CellFlattener<box_t> cell_flattener{domain_box};
};



} // namespace PHARE::core::detail


#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_GPU_HPP*/
