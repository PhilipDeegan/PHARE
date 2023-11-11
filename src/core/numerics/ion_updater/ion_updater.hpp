#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/data/ions/ions.hpp"

#include "initializer/data_provider.hpp"

#include "core/logger.hpp"

#include <cstddef>
#include <memory>

namespace PHARE::core
{
enum class UpdaterMode { domain_only = 1, all = 2 };


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField          = typename Ions::vecfield_type;
    using ParticleArray_t   = typename Ions::particle_array_type;
    using Particle_t        = typename ParticleArray_t::Particle_t;
    using PartIterator      = typename ParticleArray_t::iterator;
    using ParticleRange     = IndexRange<ParticleArray_t>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::BorisPusher<dimension, ParticleRange, Electromag, Interpolator,
                                            BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleRange, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_ = std::make_unique<Pusher>();
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
    //: pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }


    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions, GridLayout const& layout);


private:
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);



    auto interop_end(std::size_t const& it, ParticleArray_t const& particles) const
    {
        return particles.size() < EB_interop::size() ? particles.size() : EB_interop::size() * it;
    }
    auto interop_begin(std::int64_t const& end, ParticleArray_t const& particles) const
    {
        return end - EB_interop::size() >= 0 ? end - EB_interop::size() : 0;
    }
};




template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 GridLayout const& layout,
                                                                 double dt, UpdaterMode mode)
{
    PHARE_LOG_SCOPE("IonUpdater::updatePopulations");

    resetMoments(ions);
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
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions, GridLayout const& layout)
{
    ions.computeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                       Electromag const& em,
                                                                       GridLayout const& layout)
{
    PHARE_LOG_SCOPE("IonUpdater::updateAndDepositDomain_");
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    auto domainBox                = layout.AMRBox();
    auto inDomainBox              = [&](auto& particleRange) {
        return particleRange.array().partition(
            [&](auto const& cell) { return core::isIn(Point{cell}, domainBox); }, particleRange);
    };

    auto ghostBox   = grow(domainBox, partGhostWidth);
    auto inGhostBox = [&](auto& particleRange) {
        return particleRange.array().partition(
            [&](auto const& cell) { return isIn(Point{cell}, ghostBox); }, particleRange);
    };

    for (auto& pop : ions)
    {
        ParticleArray_t& domain = pop.domainParticles();

        // first push all domain particles
        // push them while still inDomainBox
        // accumulate those inDomainBox
        // erase those which left

        std::size_t deleted    = 0;
        std::size_t iterations = ceil(domain.size(), EB_interop::size());
        for (std::size_t it = iterations; it > 0; --it)
        {
            std::size_t end   = interop_end(it, domain) - deleted;
            std::size_t begin = interop_begin(end, domain);
            auto inRange      = makeIndexRange(domain, begin, end);
            auto outRange     = makeIndexRange(domain, begin, end);
            outRange          = pusher_->move_first(inRange, outRange);
            auto inDomain     = inDomainBox(
                    pusher_->move_second(outRange, em, pop.mass(), interpolator_, layout));
            interpolator_(inDomain, pop.density(), pop.flux(), layout);

            // TODO : we can erase here because we know we are working on a state
            // that has been saved in the solverPPC
            // this makes the updater quite coupled to how the solverPPC works while
            // it kind of pretends not to be by being independent object in core...
            // note we need to erase here if using the back_inserter for ghost copy
            // otherwise they will be added after leaving domain particles.

            auto toErase = makeRange(domain, inDomain.iend(), end);
            domain.erase(toErase);
            deleted += toErase.size();
        }


        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto& inputArray, bool copyInDomain = false) {
            auto outputArray       = inputArray;
            std::size_t iterations = ceil(outputArray.size(), EB_interop::size());
            for (std::size_t it = iterations; it > 0; --it)
            {
                std::size_t end      = interop_end(it, outputArray);
                std::size_t begin    = interop_begin(end, outputArray);
                auto inRange         = makeIndexRange(outputArray, begin, end);
                auto outRange        = makeIndexRange(outputArray, begin, end);
                outRange             = inGhostBox(pusher_->move_first(inRange, outRange));
                auto enteredInDomain = inDomainBox(
                    pusher_->move_second(outRange, em, pop.mass(), interpolator_, layout));
                interpolator_(enteredInDomain, pop.density(), pop.flux(), layout);
                if (copyInDomain)
                {
                    std::copy(enteredInDomain.begin(), enteredInDomain.end(),
                              std::back_inserter(domain));
                }
            }
        };

        // After this function is done domain particles overlaping ghost layers of neighbor patches
        // are sent to these neighbor's patchghost particle array.
        // After being pushed, some patch ghost particles may enter the domain. These need to be
        // copied into the domain array so they are transfered to the neighbor patch
        // ghost array and contribute to moments there too.
        // On the contrary level ghost particles entering the domain here do not need to be copied
        // since they contribute to nodes that are not shared with neighbor patches an since
        // level border nodes will receive contributions from levelghost old and new particles
        pushAndAccumulateGhosts(pop.patchGhostParticles(), true);
        pushAndAccumulateGhosts(pop.levelGhostParticles());
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
                                                                    Electromag const& em,
                                                                    GridLayout const& layout)
{
    PHARE_LOG_SCOPE("IonUpdater::updateAndDepositAll_");
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    auto domainBox                = layout.AMRBox();
    auto inDomainBox              = [&domainBox](auto& particleRange) {
        return particleRange.array().partition(
            [&](auto const& cell) { return isIn(Point{cell}, domainBox); }, particleRange);
    };

    auto ghostBox   = grow(domainBox, partGhostWidth);
    auto inGhostBox = [&](auto& particleRange) {
        return particleRange.array().partition(
            [&](auto const& cell) { return isIn(Point{cell}, ghostBox); }, particleRange);
    };

    auto inGhostLayer = [&](auto& particleRange) {
        return particleRange.array().partition(
            [&](auto const& cell) {
                return isIn(Point{cell}, ghostBox) and !isIn(Point{cell}, domainBox);
            },
            particleRange);
    };

    // push domain particles, erase from array those leaving domain
    // push patch and level ghost particles that are in ghost area (==ghost box without domain)
    // copy patch and ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in domain are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles  = pop.domainParticles();
        std::size_t deleted    = 0;
        std::size_t iterations = ceil(domainParticles.size(), EB_interop::size());
        for (std::size_t it = iterations; it > 0; --it)
        {
            std::size_t end      = interop_end(it, domainParticles) - deleted;
            std::size_t begin    = interop_begin(end, domainParticles);
            auto domainPartRange = makeIndexRange(domainParticles, begin, end);
            domainPartRange      = pusher_->move_first(domainPartRange, domainPartRange);
            auto inDomain        = inDomainBox(
                       pusher_->move_second(domainPartRange, em, pop.mass(), interpolator_, layout));
            auto toErase = makeRange(domainParticles, inDomain.iend(), end);
            domainParticles.erase(toErase);
            deleted += toErase.size();
        }


        auto pushAndCopyInDomain = [&](auto& particleArray) {
            std::size_t deleted    = 0;
            std::size_t iterations = ceil(particleArray.size(), EB_interop::size());
            for (std::size_t it = iterations; it > 0; --it)
            {
                std::size_t end    = interop_end(it, particleArray) - deleted;
                std::size_t begin  = interop_begin(end, particleArray);
                auto particleRange = makeIndexRange(particleArray, begin, end);
                particleRange      = inGhostBox(pusher_->move_first(particleRange, particleRange));
                auto inGhostLayerRange
                    = pusher_->move_second(particleRange, em, pop.mass(), interpolator_, layout);
                particleArray.export_particles(domainParticles, [&](auto const& cell) {
                    return isIn(Point{cell}, domainBox);
                });
                auto toErase = makeRange(particleArray, inGhostLayerRange.iend(), end);
                particleArray.erase(toErase);
                deleted += toErase.size();
            }
        };

        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        interpolator_(makeIndexRange(domainParticles), pop.density(), pop.flux(), layout);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
