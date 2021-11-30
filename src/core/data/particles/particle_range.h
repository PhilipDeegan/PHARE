
#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_RANGE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_RANGE_H


#include <memory>

#include "core/utilities/range/range.h"


namespace PHARE::core
{
template<typename Iterator, typename StateView>
struct ParticleRange : public Range<Iterator>
{
    using Super    = Range<Iterator>;
    using iterator = Iterator;

    ParticleRange(Range<Iterator>&& range, std::size_t pop_idx_, std::shared_ptr<StateView> view_)
        : Super{std::forward<Range<Iterator>>(range)}
        , pop_idx{pop_idx_}
        , view{view_}
    {
    }
    ParticleRange(Range<Iterator> const& range, std::size_t pop_idx_,
                  std::shared_ptr<StateView> view_)
        : Super{range}
        , pop_idx{pop_idx_}
        , view{view_}
    {
    }

    std::size_t pop_idx = -1;
    std::shared_ptr<StateView> view; // backwards relationship for particle -> mesh    };
};

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_RANGE_H */
