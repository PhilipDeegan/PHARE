#ifndef PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_CPU_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_CPU_HPP


#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <stdexcept>


namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleArrayPartitioner<AllocatorMode::CPU, /*impl = */ 0, ParticleArray>
{
public:
    using box_t    = Box<int, ParticleArray::dimension>;
    using iterator = decltype(std::declval<ParticleArray>().begin());
    using range_t  = BoxRange<box_t, iterator>;

    auto operator()(box_t const& box)
    {
        using enum LayoutMode;
        if constexpr (ParticleArray::layout_mode == SoA)
        {
            PHARE_WITH_THRUST_ELSE_THROW( //
                auto it    = SoAIteratorAdaptor::make(array);
                auto pivot = thrust::partition(
                    thrust::host, it + start, it + end,
                    [=](auto& z) { return isIn(SoAIteratorAdaptor::iCell(z), box); });
                return range_t{box, array.begin() + start,
                               array.begin() + std::distance(it, pivot)}; //
            )
        }
        // else if constexpr (any_in(ParticleArray::layout_mode, SoAVX, SoAVXTS))
        //     throw std::runtime_error("no impl");
        else if constexpr (any_in(ParticleArray::layout_mode, AoSMapped))
        {
            return array.partition(
                [&](auto& part) { return isIn(cellAsPoint(part.iCell()), box); });
        }
        else
            return partitionner(array.begin() + start, array.begin() + end, box);
    }

    // try dedupe
    auto notIn(box_t const& box)
    {
        using enum LayoutMode;
        if constexpr (ParticleArray::layout_mode == LayoutMode::SoA)
        {
            PHARE_WITH_THRUST_ELSE_THROW( //
                auto it    = SoAIteratorAdaptor::make(array);
                auto pivot = thrust::partition(
                    thrust::host, it + start, it + end,
                    [=](auto& z) { return !isIn(SoAIteratorAdaptor::iCell(z), box); });
                return range_t{box, array.begin() + start,
                               array.begin() + std::distance(it, pivot)}; //
            )
        }
        // else if constexpr (any_in(ParticleArray::layout_mode, SoAVX, SoAVXTS))
        //     throw std::runtime_error("no impl");
        else if constexpr (any_in(ParticleArray::layout_mode, AoSMapped))
        {
            return array.partition(
                [&](auto& part) { return !isIn(cellAsPoint(part.iCell()), box); });
        }
        else
            return partitionner(array.begin() + start, array.begin() + end, box,
                                [=](auto& part) { return !isIn(cellAsPoint(part.iCell()), box); });
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


template<typename ParticleArray>
class ParticleArrayPartitioner<AllocatorMode::CPU, /*impl = */ 1, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    auto operator()(box_t const& box)
    {
        static_assert(!ParticleArray::is_mapped);

        return partitionner(array.begin() + start, array.begin() + end, box,
                            [&](auto const& part) { return part.norcell >= 0; });
    }

    ParticleArray& array;
    std::size_t start = 0, end = array.size();
};



} // namespace PHARE::core::detail

#endif /*PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_CPU_HPP*/
