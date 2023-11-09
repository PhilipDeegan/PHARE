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
#include "core/data/particles/particle_array_sorting.hpp"

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

    using Box             = PHARE::core::Box<int, dimension>;
    using Interpolator    = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField        = typename Ions::vecfield_type;
    using ParticleArray_t = typename Ions::particle_array_type;
    using Particle_t      = typename ParticleArray_t::Particle_t;
    using PartIterator    = typename ParticleArray_t::const_iterator;
    using PartSorter      = ParticleCountSorting<ParticleArray_t>;

    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::BorisPusher<dimension, ParticleArray_t, Electromag, Interpolator,
                                            BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleArray_t, Electromag,
                                                 Interpolator, BoundaryCondition, GridLayout>;

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
    auto& tmpParticlesFrom(ParticleArray_t const& particles)
    {
        tmp_.update_from(particles);
        return tmp_;
    }
    auto& ghostTmpParticlesFrom(ParticleArray_t const& particles)
    {
        ghost_tmp_.update_from(particles);
        return ghost_tmp_;
    }

    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);

    std::vector<Span<Particle_t const>> ranges_;
    ParticleArray_t tmp_{Box{}}, ghost_tmp_{Box{}};
    CountingSort<ParticleArray_t, dimension> counting_sort;
    // auto ghostish_box = grow(this->particles.box());
    // // counting_sort.setup(nbrParticles, ghostish_box);
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
    auto ghostBox                 = grow(domainBox, partGhostWidth);
    auto no_op                    = [](auto& particleRange) { return particleRange; };

    for (auto& pop : ions)
    {
        ParticleArray_t& domain = pop.domainParticles();
        DEBUG_ABORT_IF(domain.size() == 0);

        pusher_->move(domain, domain, em, pop.mass(), interpolator_, layout, no_op, no_op);

        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inDomainBox


        auto pushAndAccumulateGhosts = [&](auto& inputArray, bool copyInDomain = false) {
            tmp_.clear();
            tmp_.set_as(inputArray);
            auto& outputArray = ghostTmpParticlesFrom(inputArray);

            // bleh
            pusher_->move_first(inputArray, outputArray);
            PartSorter{outputArray, counting_sort}();
            ParticleCountRangeFinder<ParticleArray_t>{outputArray}.copy_into(ghostBox, tmp_);

            pusher_->move_second(tmp_, em, pop.mass(), interpolator_, layout);
            ParticleCountSorting<ParticleArray_t>{tmp_, counting_sort}();
            ParticleCountRangeFinder<ParticleArray_t> finder{tmp_};

            auto ranges = finder.ranges(domainBox);
            if (copyInDomain) // gets interpolated on mesh later
            {
                for (auto const& range : ranges)
                    domain.insert(domain.end(), range.begin(), range.end());
            }
            else
                for (auto const& range : ranges)
                    interpolator_(range, pop.density(), pop.flux(), layout);
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


        // TODO : we can erase here because we know we are working on a state
        // that has been saved in the solverPPC
        // this makes the updater quite coupled to how the solverPPC works while
        // it kind of pretends not to be by being independent object in core...
        // note we need to erase here if using the back_inserter for ghost copy
        // otherwise they will be added after leaving domain particles.

        ParticleCountSorting<ParticleArray_t> sorter{domain, counting_sort};
        sorter();
        ParticleRangeEraser<ParticleArray_t> eraser{sorter};
        eraser.erase_outside(domainBox);
        sorter(); // possible redundant

        interpolator_(domain, pop.density(), pop.flux(), layout);

        DEBUG_ABORT_IF(domain.size() == 0);
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
    auto ghostBox                 = grow(domainBox, partGhostWidth);
    auto no_op                    = [](auto& particleRange) { return particleRange; };

    // push domain particles, erase from array those leaving domain
    // push patch and level ghost particles that are in ghost area (==ghost box without domain)
    // copy patch and ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in domain are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        ABORT_IF(domainParticles.size() == 0);

        pusher_->move(domainParticles, domainParticles, em, pop.mass(), interpolator_, layout,
                      no_op, no_op);

        auto pushAndCopyInDomain = [&](auto& particles) {
            PHARE_LOG_LINE_STR(ghostBox << " " << particles.size());
            tmp_.update_from(particles);
            particles.clear();

            pusher_->move_first(tmp_, tmp_);
            PartSorter{tmp_, counting_sort}();
            ParticleCountRangeFinder<ParticleArray_t>{tmp_}.copy_into(ghostBox, particles);

            pusher_->move_second(particles, em, pop.mass(), interpolator_, layout);

            ParticleCountSorting<ParticleArray_t> sorter{particles, counting_sort};
            sorter(); // for finding in domain and initial erase

            ParticleCountRangeFinder<ParticleArray_t> finder{particles};
            for (auto const& range : finder.ranges(domainBox))
            {
                domainParticles.insert(domainParticles.end(), range.begin(), range.end());
            }

            ParticleRangeEraser<ParticleArray_t> eraser{sorter};
            eraser.erase_outside(ghostBox);
            sorter(); // for next erase
            eraser.erase(domainBox);
            sorter(); // maybe redundant?
        };

        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        ParticleCountSorting<ParticleArray_t> sorter{domainParticles, counting_sort};
        sorter(); // sort for erase
        ParticleRangeEraser<ParticleArray_t> eraser{sorter};
        eraser.erase_outside(domainBox);

        sorter(); // for copy/stream selection

        interpolator_(domainParticles, pop.density(), pop.flux(), layout);

        DEBUG_ABORT_IF(domainParticles.size() == 0);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
