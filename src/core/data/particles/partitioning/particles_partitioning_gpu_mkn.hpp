#ifndef PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_GPU_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_GPU_HPP

#include "core/utilities/box/box.hpp"

#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>


namespace PHARE::core::detail
{

template<typename ParticleArray>
class ParticleArrayPartitioner<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>
{
public:
    using box_t    = Box<int, ParticleArray::dimension>;
    using iterator = decltype(std::declval<ParticleArray>().begin());
    using range_t  = BoxRange<box_t, iterator>;

    auto operator()(box_t const& box)
    {
        if constexpr (ParticleArray::layout_mode == LayoutMode::AoS)
        {
            auto view  = array.view();
            auto pivot = thrust::partition(
                thrust::device, view.begin() + start, view.begin() + end,
                [=] _PHARE_ALL_FN_(auto& part) { return isIn(cellAsPoint(part.iCell()), box); });
            return range_t{box, array.begin() + start,
                           array.begin() + std::distance(view.begin(), pivot)};
        }
        else if constexpr (ParticleArray::layout_mode == LayoutMode::SoA)
        {
            auto it    = SoAIteratorAdaptor::make(array);
            auto pivot = thrust::partition(
                thrust::device, it + start, it + end,
                [=] _PHARE_ALL_FN_(auto& z) { return isIn(SoAIteratorAdaptor::iCell(z), box); });
            return range_t{box, array.begin() + start, array.begin() + std::distance(it, pivot)};
        }
    }


    template<std::size_t S>
    auto operator()(std::array<box_t, S> const& boxes)
    {
        static_assert(S > 0);
        auto iterators = generate_from(
            [&](auto const& box) { return range_t{box, array.begin(), array.begin()}; }, boxes);
        for (std::size_t i = 0; i < boxes.size(); ++i)
            start += (iterators[i] = (*this)(boxes[i])).size();
        return iterators;
    }

    auto operator()(std::vector<box_t> const& boxes)
    {
        std::vector<range_t> iterators;
        for (auto const& box : boxes)
            start += iterators.emplace_back((*this)(box)).size();
        return iterators;
    }

    ParticleArray& array;
    std::size_t start = 0, end = array.size();
};



} // namespace PHARE::core::detail


#endif /*PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_GPU_HPP*/
