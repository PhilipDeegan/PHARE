#ifndef PHARE_CORE_LLNL_ION_UPDATER_H
#define PHARE_CORE_LLNL_ION_UPDATER_H

#include "core/numerics/ion_updater/ion_updater.h"
// TODO alpha coef for interpolating new and old levelGhost should be given somehow...


namespace PHARE::core::llnl
{
template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box          = PHARE::core::Box<int, dimension>;
    using Interpolator = PHARE::core::Interpolator<dimension, interp_order>;

    using VecField          = typename Ions::vecfield_type;
    using ParticleArray     = typename Ions::particle_array_type;
    using PartIterator      = typename ParticleArray::iterator;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher            = PHARE::core::Pusher<dimension, PartIterator, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }

    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::all);

    void updateIons(Ions& ions, GridLayout const& layout);

private:
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);
    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_{};
};




template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout, Offloader>::updatePopulations(
    Ions& ions, Electromag const& em, GridLayout const& layout, double dt, UpdaterMode mode)
{
    pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
    {
        updateAndDepositDomain_(ions, em, layout);
    }
    else
    {
        updateAndDepositAll_(ions, em, layout);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout, Offloader>::updateIons(Ions& ions,
                                                                     GridLayout const& layout)
{
    fixMomentGhosts(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout, Offloader>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonUpdater<Ions, Electromag, GridLayout, Offloader>::updateAndDepositDomain_(
    Ions& ions, Electromag const& em, GridLayout const& layout)
{
    // domain particles are offloaded so are not pushed/interpoloated here.

    auto domainBox = layout.AMRBox();

    auto inDomainBox = [&domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, domainBox);
    };

    auto constexpr partGhostWidth = GridLayout::ghostWidthForParticles();
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

    auto ghostSelector = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };

    for (auto& pop : ions)
    {
        ParticleArray& domain = pop.domainParticles();
        ParticleArray tmpPatchGhost;
        ParticleArray tmpLevelGhost;

        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // some will leave the ghost area
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto& inputArray, auto& outputArray,
                                           bool copyInDomain = false) {
            outputArray.resize(inputArray.size());

            auto inRange  = makeRange(std::begin(inputArray), std::end(inputArray));
            auto outRange = makeRange(std::begin(outputArray), std::end(outputArray));

            auto firstGhostOut = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_,
                                               ghostSelector, layout);

            auto endInDomain = std::partition(firstGhostOut, std::end(outputArray), inDomainBox);

            interpolator_(firstGhostOut, endInDomain, pop.density(), pop.flux(), layout);

            if (copyInDomain)
                std::copy(firstGhostOut, endInDomain, std::back_inserter(domain));
        };

        // After this function is done domain particles overlaping ghost layers of neighbor patches
        // are sent to these neighbor's patchghost particle array.
        // After being pushed, some patch ghost particles may enter the domain. These need to be
        // copied into the domain array so they are transfered to the neighbor patch
        // ghost array and contribute to moments there too.
        // On the contrary level ghost particles entering the domain here do not need to be copied
        // since they contribute to nodes that are not shared with neighbor patches an since
        // level border nodes will receive contributions from levelghost old and new particles
        pushAndAccumulateGhosts(pop.patchGhostParticles(), tmpPatchGhost, true);
        pushAndAccumulateGhosts(pop.levelGhostParticles(), tmpLevelGhost);
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout, Offloader>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonUpdater<Ions, Electromag, GridLayout, Offloader>::updateAndDepositAll_(
    Ions& ions, Electromag const& em, GridLayout const& layout)
{
    auto constexpr partGhostWidth = GridLayout::ghostWidthForParticles();
    auto domainBox                = layout.AMRBox();
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

    auto ghostSelector = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };

    auto inDomainSelector
        = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };

    // push domain particles, erase from array those leaving domain
    // push patch and level ghost particles that are in ghost area (==ghost box without domain)
    // copy patch and ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in domain are to be interpolated on mesh.

    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto firstOutside     = domainParticles.end(); // TODO FIX
        domainParticles.erase(firstOutside, std::end(domainParticles));
        assert(domainParticles.size());

        auto pushAndCopyInDomain = [&](auto& particleArray) {
            auto range = makeRange(particleArray);

            auto firstOutGhostBox
                = pusher_->move(range, em, pop.mass(), interpolator_, ghostSelector, layout);

            std::copy_if(firstOutGhostBox, std::end(particleArray),
                         std::back_inserter(domainParticles), inDomainSelector);

            particleArray.erase(firstOutGhostBox, std::end(particleArray));
        };

        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        // only interpolate new domain particles here
        interpolator_(firstOutside, std::end(domainParticles), pop.density(), pop.flux(), layout);
    }
}



} // namespace PHARE::core::llnl


#endif // ION_UPDATER_H
