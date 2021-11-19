#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H


#include "core/utilities/box/box.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/moments/moments.h"

#include "core/data/ions/ions.h"

#include "initializer/data_provider.h"

#include "core/logger.h"

#include <memory>

// TODO alpha coef for interpolating new and old levelGhost should be given somehow...


namespace PHARE::core
{
enum class UpdaterMode { domain_only = 1, all = 2 };

template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr std::size_t def_op_size = 1e6;

    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
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
    Interpolator interpolator_;

    std::size_t operating_particle_size(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict.contains("operating_particle_size"))
            return dict["operating_particle_size"].template to<std::size_t>();
        return def_op_size;
    }

public:
    IonUpdater(std::string pusher_name, std::size_t operating_particle_size = def_op_size)
        : pusher_{makePusher(pusher_name)}
        , particle_EBs(operating_particle_size)
    {
    }

    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : IonUpdater{dict["pusher"]["name"].template to<std::string>(),
                     operating_particle_size(dict)}
    {
    }

    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions, GridLayout const& layout);


private:
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);


    // std::vector<tuple_fixed_type<
    using EB = tuple_fixed_type<double, 3>;
    std::vector<tuple_fixed_type<EB, 2>> particle_EBs;

    std::size_t partition_idx = 0, domain_particles = 0;
    std::vector<std::pair<std::size_t, std::size_t>> partition_from_to;
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
    fixMomentGhosts(ions, layout);
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

    auto domainBox = layout.AMRBox();

    auto inDomainBox
        = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };

    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

    auto inGhostLayer = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };
    auto inGhostBox = [&ghostBox](auto& part) { return core::isIn(cellAsPoint(part), ghostBox); };

    for (auto& pop : ions)
    {
        ParticleArray& domain = pop.domainParticles();

        // first push all domain particles
        // push them while still inDomainBox
        // accumulate those inDomainBox

        auto inRange  = makeRange(domain);
        auto outRange = makeRange(domain);

        auto newEnd = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_, inDomainBox,
                                    layout, particle_EBs);

        interpolator_(std::begin(domain), newEnd, pop.density(), pop.flux(), layout);

        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // some will leave the ghost area
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto& inputArray, bool copyInDomain = false) {
            ParticleArray outputArray(inputArray.size());

            auto inRange  = makeRange(inputArray);
            auto outRange = makeRange(outputArray);

            auto firstGhostOut = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_,
                                               inGhostBox, inGhostLayer, layout, particle_EBs);

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
        // KLOG(INF) << "move1 patchGhostParticles";
        pushAndAccumulateGhosts(pop.patchGhostParticles(), true);
        // KLOG(INF) << "move1 levelGhostParticles";
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
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

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

    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();

        auto domainPartRange = makeRange(domainParticles);

        auto firstOutside = pusher_->move(domainPartRange, domainPartRange, em, pop.mass(),
                                          interpolator_, inDomainBox, layout, particle_EBs);

        domainParticles.erase(firstOutside, std::end(domainParticles));

        auto pushAndCopyInDomain = [&](auto& particleArray) {
            auto range            = makeRange(particleArray);
            auto firstOutGhostBox = pusher_->move(range, range, em, pop.mass(), interpolator_,
                                                  inGhostBox, inGhostLayer, layout, particle_EBs);

            std::copy_if(firstOutGhostBox, std::end(particleArray),
                         std::back_inserter(domainParticles), inDomainBox);

            particleArray.erase(firstOutGhostBox, std::end(particleArray));
        };


        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        interpolator_(std::begin(domainParticles), std::end(domainParticles), pop.density(),
                      pop.flux(), layout);
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_H




// #ifndef PHARE_ION_UPDATER_H
// #define PHARE_ION_UPDATER_H


// #include "core/utilities/box/box.h"
// #include "core/numerics/interpolator/interpolator.h"
// #include "core/numerics/pusher/pusher.h"
// #include "core/numerics/pusher/pusher_factory.h"
// #include "core/numerics/boundary_condition/boundary_condition.h"
// #include "core/numerics/moments/moments.h"

// #include "core/data/ions/ions.h"

// #include "initializer/data_provider.h"

// #include "core/logger.h"

// #include <memory>

// // TODO alpha coef for interpolating new and old levelGhost should be given somehow...


// namespace PHARE::core
// {
// enum class UpdaterMode { domain_only = 1, all = 2 };

// template<typename Ions, typename Electromag, typename GridLayout>
// class IonUpdater
// {
// public:
//     static constexpr std::size_t def_op_size = 1e6;

//     static constexpr auto dimension    = GridLayout::dimension;
//     static constexpr auto interp_order = GridLayout::interp_order;

//     using Box               = PHARE::core::Box<int, dimension>;
//     using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
//     using VecField          = typename Ions::vecfield_type;
//     using ParticleArray     = typename Ions::particle_array_type;
//     using PartIterator      = typename ParticleArray::iterator;
//     using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
//     using Pusher            = PHARE::core::Pusher<dimension, PartIterator, Electromag,
//     Interpolator,
//                                        BoundaryCondition, GridLayout>;

// private:
//     constexpr static auto makePusher
//         = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag,
//         Interpolator,
//                                                  BoundaryCondition, GridLayout>;

//     std::unique_ptr<Pusher> pusher_;
//     Interpolator interpolator_;

//     std::size_t operating_particle_size(PHARE::initializer::PHAREDict const& dict)
//     {
//         if (dict.contains("operating_particle_size"))
//             return dict["operating_particle_size"].template to<std::size_t>();
//         return def_op_size;
//     }

// public:
//     IonUpdater(std::string pusher_name, std::size_t operating_particle_size = def_op_size)
//         : pusher_{makePusher(pusher_name)}
//         , particle_EBs(operating_particle_size)
//     {
//     }

//     IonUpdater(PHARE::initializer::PHAREDict const& dict)
//         : IonUpdater{dict["pusher"]["name"].template to<std::string>(),
//                      operating_particle_size(dict)}
//     {
//     }

//     void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
//                            UpdaterMode = UpdaterMode::all);


//     void updateIons(Ions& ions, GridLayout const& layout);


// private:
//     void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);

//     void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);

//     template<typename Population, typename Selector>
//     void push_domain(Population& pop, Selector& selector, Electromag const& em,
//                      GridLayout const& layout);

//     template<typename Population, typename GhostParticles, typename GhostParticleRange>
//     void copy_to_domain(Population& pop, GhostParticles&, GhostParticleRange in_domain,
//                         bool clean_domain = false);

//     template<typename Particles>
//     void _clean_domain(Particles& domain);

//     // std::vector<tuple_fixed_type<
//     using EB = tuple_fixed_type<double, 3>;
//     std::vector<tuple_fixed_type<EB, 2>> particle_EBs;

//     std::size_t partition_idx = 0, domain_particles = 0;
//     std::vector<std::pair<std::size_t, std::size_t>> partition_from_to;
// };




// template<typename Ions, typename Electromag, typename GridLayout>
// void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const&
// em,
//                                                                  GridLayout const& layout,
//                                                                  double dt, UpdaterMode mode)
// {
//     PHARE_LOG_SCOPE("IonUpdater::updatePopulations");

//     resetMoments(ions);
//     pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

//     if (mode == UpdaterMode::domain_only)
//     {
//         updateAndDepositDomain_(ions, em, layout);
//     }
//     else
//     {
//         updateAndDepositAll_(ions, em, layout);
//     }
// }



// template<typename Ions, typename Electromag, typename GridLayout>
// void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions, GridLayout const& layout)
// {
//     fixMomentGhosts(ions, layout);
//     ions.computeDensity();
//     ions.computeBulkVelocity();
// }


// template<typename Ions, typename Electromag, typename GridLayout>
// template<typename Population, typename Selector>
// void IonUpdater<Ions, Electromag, GridLayout>::push_domain(Population& pop, Selector& selector,
//                                                            Electromag const& em,
//                                                            GridLayout const& layout)
// {
//     auto& particles = pop.domainParticles();

//     std::size_t remaining = particles.size();
//     auto curr_iter        = std::begin(particles);
//     auto& ebs             = particle_EBs;

//     {
//         std::size_t needs = 1;
//         if (remaining > ebs.size())
//             needs = remaining / ebs.size() + 1;
//         if (partition_from_to.size() < needs)
//             partition_from_to.resize(needs);
//     }

//     partition_idx = 0;
//     while (remaining > 0)
//     {
//         std::size_t operate = remaining >= ebs.size() ? ebs.size() : remaining;
//         auto range          = makeRange(curr_iter, curr_iter + operate);
//         auto newEnd
//             = pusher_->move(range, range, em, pop.mass(), interpolator_, selector, layout, ebs);
//         interpolator_(curr_iter, newEnd, pop.density(), pop.flux(), layout);
//         curr_iter += operate;
//         remaining -= operate;

//         auto particles_remaining = std::distance(range.begin(), newEnd);
//         auto relative_begin      = ebs.size() * partition_idx;

//         partition_from_to[partition_idx]
//             = {relative_begin + particles_remaining, operate - particles_remaining};

//         // KLOG(INF) << partition_from_to[partition_idx].first << " "        <<
//         // partition_from_to[partition_idx].second;

//         ++partition_idx;
//     }
//     partition_idx = 0;
// }


// template<typename Ions, typename Electromag, typename GridLayout>
// template<typename Population, typename GhostParticles, typename GhostParticleRange>
// void IonUpdater<Ions, Electromag, GridLayout>::copy_to_domain(Population& pop,
//                                                               GhostParticles& ghosts,
//                                                               GhostParticleRange new_in_domain,
//                                                               bool clean_domain_)
// {
//     auto& domain                     = pop.domainParticles();
//     auto const& copy_from            = partition_from_to[0].first;
//     auto const& n_domain_replaceable = partition_from_to[0].second;
//     auto n_new_domain_particles      = new_in_domain.size();

//     // KLOG(INF) << n_new_domain_particles;
//     // KLOG(INF) << n_domain_replaceable;

//     if (n_new_domain_particles > 0)
//     {
//         if (n_domain_replaceable >= n_new_domain_particles)
//         {
//             // KLOG(INF);
//             std::copy(new_in_domain.begin(), new_in_domain.end(), domain.begin() + copy_from);

//             partition_from_to[0].first += n_new_domain_particles;
//             partition_from_to[0].second -= n_new_domain_particles;
//         }
//         else
//         {
//             // KLOG(INF);
//             // replace replaceable before back insertion
//             if (n_domain_replaceable > 0)
//                 std::copy(new_in_domain.begin(), new_in_domain.begin() + n_domain_replaceable,
//                           domain.begin() + copy_from);

//             std::copy(new_in_domain.begin() + n_domain_replaceable, new_in_domain.end(),
//                       std::back_inserter(domain));

//             partition_from_to[0].first += partition_from_to[0].second;
//             partition_from_to[0].second = 0; // kill
//         }
//     }

//     // KLOG(INF) << domain.size();
//     if (clean_domain_)
//         _clean_domain(domain);
//     // KLOG(INF) << domain.size();
// }

// template<typename Ions, typename Electromag, typename GridLayout>
// template<typename Particles>
// void IonUpdater<Ions, Electromag, GridLayout>::_clean_domain(Particles& domain)
// {
//     // KLOG(INF);

//     auto const& n_domain_deletable = partition_from_to[0].second;

//     if (n_domain_deletable > 0)
//     {
//         auto const& delete_from = partition_from_to[0].first;
//         auto const delete_to    = delete_from + n_domain_deletable;
//         domain.erase(domain.begin() + delete_from, domain.begin() + delete_to);
//     }

//     // KLOG(INF) << domain.size();
// }


// template<typename Ions, typename Electromag, typename GridLayout>
// /**
//  * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
//    evolves moments from time n to n+1 without updating particles, which stay at time n
//  */
// void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
//                                                                        Electromag const& em,
//                                                                        GridLayout const& layout)
// {
//     PHARE_LOG_SCOPE("IonUpdater::updateAndDepositDomain_");

//     auto domainBox = layout.AMRBox();

//     auto inDomainBox
//         = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };

//     auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
//     auto ghostBox{domainBox};
//     ghostBox.grow(partGhostWidth);

//     auto inGhostLayer = [&ghostBox, &domainBox](auto const& part) {
//         auto cell = cellAsPoint(part);
//         return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
//     };
//     auto inGhostBox = [&ghostBox](auto& part) { return core::isIn(cellAsPoint(part), ghostBox);
//     };


//     // first push all domain particles
//     // push them while still inDomainBox
//     // accumulate those inDomainBox

//     for (auto& pop : ions)
//         push_domain(pop, inDomainBox, em, layout);


//     for (auto& pop : ions)
//     {
//         // then push patch and level ghost particles
//         // push those in the ghostArea (i.e. stop pushing if they're not out of it)
//         // some will leave the ghost area
//         // deposit moments on those which leave to go inDomainBox

//         auto pushAndAccumulateGhosts = [&](auto& inputArray, bool copyInDomain = false) {
//             ParticleArray outputArray(inputArray.size());

//             auto inRange  = makeRange(inputArray);
//             auto outRange = makeRange(outputArray);

//             auto firstGhostOut = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_,
//                                                inGhostBox, inGhostLayer, layout, particle_EBs);

//             auto endInDomain = std::partition(firstGhostOut, std::end(outputArray), inDomainBox);

//             interpolator_(firstGhostOut, endInDomain, pop.density(), pop.flux(), layout);

//             if (copyInDomain)
//                 copy_to_domain(pop, outputArray, makeRange(firstGhostOut, endInDomain));
//         };

//         // After this function is done domain particles overlaping ghost layers of neighbor
//         patches
//         // are sent to these neighbor's patchghost particle array.
//         // After being pushed, some patch ghost particles may enter the domain. These need to be
//         // copied into the domain array so they are transfered to the neighbor patch
//         // ghost array and contribute to moments there too.
//         // On the contrary level ghost particles entering the domain here do not need to be
//         copied
//         // since they contribute to nodes that are not shared with neighbor patches an since
//         // level border nodes will receive contributions from levelghost old and new particles
//         // KLOG(INF) << "move1 patchGhostParticles";
//         pushAndAccumulateGhosts(pop.patchGhostParticles(), true);
//         // KLOG(INF) << "move1 levelGhostParticles";
//         pushAndAccumulateGhosts(pop.levelGhostParticles());
//     }
// }


// template<typename Ions, typename Electromag, typename GridLayout>
// /**
//  * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
//    evolves moments and particles from time n to n+1
//  */
// void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
//                                                                     Electromag const& em,
//                                                                     GridLayout const& layout)
// {
//     PHARE_LOG_SCOPE("IonUpdater::updateAndDepositAll_");

//     auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
//     auto domainBox                = layout.AMRBox();
//     auto ghostBox{domainBox};
//     ghostBox.grow(partGhostWidth);

//     auto inGhostLayer = [&ghostBox, &domainBox](auto const& part) {
//         auto cell = cellAsPoint(part);
//         return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
//     };
//     auto inGhostBox = [&ghostBox](auto& part) { return core::isIn(cellAsPoint(part), ghostBox);
//     }; auto inDomainBox
//         = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };

//     // push domain particles, erase from array those leaving domain
//     // push patch and level ghost particles that are in ghost area (==ghost box without domain)
//     // copy patch and ghost particles out of ghost area that are in domain, in particle array
//     // finally all particles in domain are to be interpolated on mesh.

//     for (auto& pop : ions)
//         push_domain(pop, inDomainBox, em, layout);

//     for (auto& pop : ions)
//     {
//         auto& domainParticles = pop.domainParticles();

//         // domainParticles.erase(firstOutside, std::end(domainParticles));

//         auto pushAndCopyInDomain = [&](auto& ghostParticles, bool clean_domain = false) {
//             auto range         = makeRange(ghostParticles);
//             auto firstGhostOut = pusher_->move(range, range, em, pop.mass(), interpolator_,
//                                                inGhostBox, inGhostLayer, layout, particle_EBs);

//             auto endInDomain = std::partition(firstGhostOut, std::end(ghostParticles),
//             inDomainBox);

//             copy_to_domain(pop, ghostParticles, makeRange(firstGhostOut, endInDomain),
//                            clean_domain);

//             // std::copy_if(firstGhostOut, std::end(ghostParticles),
//             //              std::back_inserter(domainParticles), inDomainBox);

//             ghostParticles.erase(firstGhostOut, std::end(ghostParticles));
//         };

//         // KLOG(INF) << "move2 patchGhostParticles";
//         pushAndCopyInDomain(pop.patchGhostParticles());
//         // KLOG(INF) << "move2 levelGhostParticles";
//         pushAndCopyInDomain(pop.levelGhostParticles(), /*clean_domain=*/true);

//         interpolator_(domainParticles, pop.density(), pop.flux(), layout);
//     }
// }



// } // namespace PHARE::core


// #endif // ION_UPDATER_H
