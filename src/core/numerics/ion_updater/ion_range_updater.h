

#ifndef PHARE_ION_RANGE_UPDATER_H
#define PHARE_ION_RANGE_UPDATER_H


#include "core/utilities/box/box.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/moments/moments.h"

#include "core/data/ions/ions.h"

#include "initializer/data_provider.h"
#include "core/utilities/range/range_replacer.h"

#include "core/logger.h"

#include <thread>
#include <mutex>
#include <memory>

#include "core/numerics/ion_updater/ion_updater.h" // UpdaterMode
// TODO alpha coef for interpolating new and old levelGhost should be given somehow...


namespace PHARE::core
{
// enum class UpdaterMode { domain_only = 1, all = 2 };

template<typename HybridModel, typename GridLayout>
class IonRangeUpdater
{
public:
    static constexpr bool atomic_interp      = true;
    static constexpr std::size_t def_op_size = 1e6;

    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Ions = typename HybridModel::ions_type;

    using HybridStateView = typename HybridModel::StateView_t;
    using Electromag      = typename HybridStateView::Electromag_t;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order, atomic_interp>;
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

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_{};

    std::size_t static operating_particle_size(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict.contains("operating_particle_size"))
            return dict["operating_particle_size"].template to<std::size_t>();
        return def_op_size;
    }

public:
    IonRangeUpdater(RangeSynchrotron<ParticleArray>& synchrotron_, std::uint16_t thread_idx,
                    std::string pusher_name, std::size_t operating_particle_size = def_op_size)
        : pusher_{makePusher(pusher_name)}
        , synchrotron{synchrotron_}
        , thread_idx_{thread_idx}
        , particle_EBs(operating_particle_size)
    {
    }

    IonRangeUpdater(RangeSynchrotron<ParticleArray>& synchrotron_, std::uint16_t thread_idx,
                    PHARE::initializer::PHAREDict const& dict)
        : IonRangeUpdater{synchrotron_, thread_idx,
                          dict["pusher"]["name"].template to<std::string>(),
                          operating_particle_size(dict)}
    {
    }

    template<typename ParticleRanges>
    void updatePopulations(ParticleRanges& particleRanges, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions, GridLayout const& layout);


private:
    template<typename ParticleRanges>
    void updateAndDepositDomain_(ParticleRanges& particleRanges);

    template<typename ParticleRanges>
    void updateAndDepositAll_(ParticleRanges& particleRanges);

    template<typename ParticleRanges, typename Selector>
    void push_domain(ParticleRanges& particleRanges, Selector& selector);


    /** Used here **/
    RangeSynchrotron<ParticleArray>& synchrotron; // Only one
    std::uint16_t thread_idx_ = -1;
    std::unique_ptr<RangeReplacer<ParticleArray>> refiller; // one of many
    /** Used here **/

    /** Used between boris and interpolator **/
    using EB = tuple_fixed_type<double, 3>;
    std::vector<tuple_fixed_type<EB, 2>> particle_EBs;
    /** Used between boris and interpolator **/
};




template<typename HybridModel, typename GridLayout>
template<typename ParticleRanges>
void IonRangeUpdater<HybridModel, GridLayout>::updatePopulations(ParticleRanges& particleRanges,
                                                                 double dt, UpdaterMode mode)
{
    PHARE_LOG_SCOPE("IonRangeUpdater::updatePopulations");

    for (auto& [layout, daos] : particleRanges)
    {
        pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

        if (mode == UpdaterMode::domain_only)
        {
            for (auto& dao : daos)
                updateAndDepositDomain_(dao);
        }
        else
        {
            for (auto& dao : daos)
                updateAndDepositAll_(dao);
        }
    }
}



template<typename HybridModel, typename GridLayout>
void IonRangeUpdater<HybridModel, GridLayout>::updateIons(Ions& ions, GridLayout const& layout)
{
    fixMomentGhosts(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}


template<typename HybridModel, typename GridLayout>
template<typename Particles, typename Selector>
void IonRangeUpdater<HybridModel, GridLayout>::push_domain(Particles& particles, Selector& selector)
{
    auto& layout        = particles.domain.view->layout;
    auto& em            = particles.domain.view->electromag;
    auto& pop           = *(particles.domain.view->ions.begin() + particles.domain.pop_idx);
    auto const pop_size = (*pop.domain).size();
    refiller
        = std::make_unique<RangeReplacer<ParticleArray>>(*pop.domain, synchrotron, thread_idx_);

    for (auto& range : ranges(particles.domain, particle_EBs.size()))
    {
        auto rangeOut = range; // is modified in boris
        auto newEnd = pusher_->move(range, rangeOut, em, pop.mass, interpolator_, selector, layout,
                                    particle_EBs);

        interpolator_(range.begin(), newEnd, pop.density, pop.flux, layout);

        refiller->add_replaceable(range, newEnd);
    }
}




template<typename HybridModel, typename GridLayout>
template<typename Particles>
/**
 * @brief IonRangeUpdater<HybridModel, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonRangeUpdater<HybridModel, GridLayout>::updateAndDepositDomain_(Particles& particles)
{
    PHARE_LOG_SCOPE("IonRangeUpdater::updateAndDepositDomain_");

    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    constexpr bool copy_src       = true; // ghost ParticleArray outputArray is temporary

    auto& layout = particles.domain.view->layout;
    auto& em     = particles.domain.view->electromag;
    auto& pop    = *(particles.domain.view->ions.begin() + particles.domain.pop_idx);

    auto domainBox = layout.AMRBox();
    auto ghostBox  = domainBox.copy().grow(partGhostWidth);

    auto inDomainBox
        = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };
    auto inGhostLayer = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };
    auto inGhostBox = [&ghostBox](auto& part) { return core::isIn(cellAsPoint(part), ghostBox); };


    // first push all domain particles push them while
    // still inDomainBox accumulate those inDomainBox

    push_domain(particles, inDomainBox);

    // then push patch and level ghost particles
    // push those in the ghostArea (i.e. stop pushing if they're not out of it)
    // some will leave the ghost area
    // deposit moments on those which leave to go inDomainBox

    auto pushAndAccumulateGhosts = [&](auto& inputArrayPtr, bool copyInDomain = false) {
        if (!inputArrayPtr)
            return;

        auto& inputArray = *inputArrayPtr;
        ParticleArray outputArray(inputArray.size());

        auto inRange  = makeRange(inputArray);
        auto outRange = makeRange(outputArray);

        auto firstGhostOut = pusher_->move(inRange, outRange, em, pop.mass, interpolator_,
                                           inGhostBox, inGhostLayer, layout, particle_EBs);

        auto endInDomain = std::partition(firstGhostOut, std::end(outputArray), inDomainBox);

        interpolator_(firstGhostOut, endInDomain, pop.density, pop.flux, layout);

        if (copyInDomain)
            refiller->replace(outputArray, makeRange(firstGhostOut, endInDomain), copy_src);
    };


    // After this function is done domain particles overlaping ghost layers of neighbor patches
    // are sent to these neighbor's patchghost particle array.
    // After being pushed, some patch ghost particles may enter the domain. These need to be
    // copied into the domain array so they are transfered to the neighbor patch
    // ghost array and contribute to moments there too.
    // On the contrary level ghost particles entering the domain here do not need to be copied
    // since they contribute to nodes that are not shared with neighbor patches an since
    // level border nodes will receive contributions from levelghost old and new
    // particles


    pushAndAccumulateGhosts(particles.level_ghost);
    pushAndAccumulateGhosts(particles.patch_ghost, /*copyInDomain=*/true);

    refiller->erase();
}


template<typename HybridModel, typename GridLayout>
template<typename Particles>
/**
 * @brief IonRangeUpdater<HybridModel, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonRangeUpdater<HybridModel, GridLayout>::updateAndDepositAll_(Particles& particles)
{
    PHARE_LOG_SCOPE("IonRangeUpdater::updateAndDepositAll_");

    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    auto constexpr copy_src       = false; // no need here

    auto& layout = particles.domain.view->layout;
    auto& em     = particles.domain.view->electromag;
    auto& pop    = *(particles.domain.view->ions.begin() + particles.domain.pop_idx);

    auto domainBox = layout.AMRBox();
    auto ghostBox  = domainBox.copy().grow(partGhostWidth);

    auto inGhostLayer = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };
    auto inGhostBox = [&ghostBox](auto& part) { return core::isIn(cellAsPoint(part), ghostBox); };
    auto inDomainBox
        = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };


    // push domain particles, erase from array those leaving domain
    // push patch and level ghost particles that are in ghost area (==ghost box without domain)
    // copy patch and ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in domain are to be interpolated on mesh.

    push_domain(particles, inDomainBox);

    auto pushAndCopyInDomain = [&](auto& ghostParticlesPtr, auto& ghost_src,
                                   bool clean_src = false) {
        if (!ghostParticlesPtr)
            return;

        auto& ghostParticles = *ghostParticlesPtr;
        auto range           = makeRange(ghostParticles);
        auto firstGhostOut   = pusher_->move(range, range, em, pop.mass, interpolator_, inGhostBox,
                                           inGhostLayer, layout, particle_EBs);

        auto endInDomain = std::partition(firstGhostOut, std::end(ghostParticles), inDomainBox);

        interpolator_(firstGhostOut, endInDomain, pop.density, pop.flux, layout);

        refiller->replace(ghost_src, makeRange(firstGhostOut, endInDomain), copy_src, clean_src);
    };

    pushAndCopyInDomain(particles.level_ghost, *pop.level_ghost);
    pushAndCopyInDomain(particles.patch_ghost, *pop.patch_ghost, /*clean_src=*/true);

    refiller->erase();
}


} // namespace PHARE::core


#endif // ION_UPDATER_H
