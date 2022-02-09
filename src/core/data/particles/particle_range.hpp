
#ifndef PHARE_CORE_DATA_PARTICLE_RANGE_HPP
#define PHARE_CORE_DATA_PARTICLE_RANGE_HPP

#include "core/utilities/range/range.hpp"


namespace PHARE::core
{
template<typename Iterator>
struct ParticleRange : public Range<Iterator>
{
    auto static constexpr is_contiguous = Iterator::is_contiguous;

    using Super    = Range<Iterator>;
    using iterator = Iterator;

    ParticleRange(Iterator begin, Iterator end)
        : Super{begin, end}
    {
    }
};


template<typename Iterator>
ParticleRange<Iterator> makeParticleRange(Iterator begin, Iterator last)
{
    return {std::forward<Iterator>(begin), std::forward<Iterator>(last)};
}

template<typename Container>
auto makeParticleRange(Container& container)
{
    return makeParticleRange(std::begin(container), std::end(container));
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLE_RANGE_HPP */
