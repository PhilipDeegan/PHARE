#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include <cstddef>
#include <functional>

namespace PHARE
{
namespace core
{
    template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
             typename BoundaryCondition, typename GridLayout, bool UseParticleIndex = true>
    class Pusher
    {
    protected:
        static auto constexpr dimension = GridLayout::dimension;

        using ParticleSelector = std::function<ParticleRange(ParticleRange&)>;

    public:
        // TODO : to really be independant on boris which has 2 push steps
        // we should have an arbitrary number of selectors, 1 per push step
        virtual ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                                   Electromag const& emFields, floater_t<4> mass,
                                   Interpolator& interpolator, GridLayout const& layout,
                                   ParticleSelector firstSelector, ParticleSelector secondSelector)
            = 0;


        virtual void setMeshAndTimeStep(std::array<floater_t<4>, dim> const ms,
                                        floater_t<4> const ts)
            = 0;

        virtual ~Pusher() {}
    };

} // namespace core
} // namespace PHARE

#endif
