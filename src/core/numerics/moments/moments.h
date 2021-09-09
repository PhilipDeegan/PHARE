#ifndef MOMENTS_H
#define MOMENTS_H

#include <iterator>

#include "core/numerics/interpolator/interpolator.h"


namespace PHARE::core
{
struct DomainDeposit
{
};

struct PatchGhostDeposit
{
};

struct LevelGhostDeposit
{
};

template<typename Ions, typename GridLayout, typename DepositTag>
void deposit_particles(Ions& ions, GridLayout& layout,
                       Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate,
                       DepositTag)
{
    for (auto& pop : ions)
    {
        auto& density = pop.density();
        auto& flux    = pop.flux();

        if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            interpolate(pop.domainParticles(), density, flux, layout);

        else if constexpr (std::is_same_v<DepositTag, PatchGhostDeposit>)
            interpolate(pop.patchGhostParticles(), density, flux, layout);

        else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            interpolate(pop.levelGhostParticlesOld(), density, flux, layout);

        else
            assert(false);
    }
}

#if defined(HAVE_UMPIRE)
template<typename Ions, typename GridLayout, typename DepositTag>
void deposit_umpire_particles(
    Ions& ions, GridLayout& layout,
    Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate, DepositTag)
{
    // assert(false);
}
#endif

} // namespace PHARE::core



namespace PHARE
{
namespace core
{
    template<typename Ions>
    void resetMoments(Ions& ions)
    {
        for (auto& pop : ions)
        {
            pop.density().zero();
            pop.flux().zero();
        }
    }


    template<typename Ions, typename GridLayout, typename DepositTag>
    void depositParticles(Ions& ions, GridLayout& layout,
                          Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate,
                          DepositTag tag)
    {
        using FieldT = typename Ions::field_type;

        if constexpr (FieldT::is_host_mem)
            deposit_particles(ions, layout, interpolate, tag);
#if defined(HAVE_UMPIRE)
        else
            deposit_umpire_particles(ions, layout, interpolate, tag);
#else
        else
            assert(false);
#endif
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_H
