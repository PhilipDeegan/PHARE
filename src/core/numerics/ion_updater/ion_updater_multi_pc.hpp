#ifndef PHARE_ION_UPDATER_PC_HPP
#define PHARE_ION_UPDATER_PC_HPP

#include "ion_updater_def.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

namespace PHARE::core
{
template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdaterMultiPC
{
    using Particles                     = typename Ions::particle_array_type;
    bool constexpr static atomic_interp = Particles::alloc_mode == AllocatorMode::GPU_UNIFIED;
    // static_assert(Particles::alloc_mode == AllocatorMode::GPU_UNIFIED);

public:
    auto static constexpr dimension    = GridLayout::dimension;
    auto static constexpr interp_order = GridLayout::interp_order;
    using Box_t                        = Box<int, dimension>;
    using Interpolator_t               = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t              = Interpolating<Particles, interp_order, atomic_interp>;
    // using Particle_t    = typename Particles::Particle_t;

    using Pusher_t = MultiBorisPusher<dimension, Interpolator_t, GridLayout>;

private:
    Interpolating_t interpolator_;


public:
    IonUpdaterMultiPC(PHARE::initializer::PHAREDict const& /*dict*/) {}

    template<typename ModelViews>
    void updatePopulations(ModelViews& views, double const& dt, UpdaterMode = UpdaterMode::all);


    static void updateIons(Ions& ions)
    {
        // fixMomentGhosts(ions, layout);
        ions.computeDensity();
        ions.computeBulkVelocity();
    }

    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(make_particles<Particles>(Box_t{}));
    }


private:
    template<typename ModelViews>
    void updateAndDepositDomain_(ModelViews& views);

    template<typename ModelViews>
    void updateAndDepositAll_(ModelViews& views);



    // dealloced on regridding/load balancing coarsest
    Particles tmp_particles_;

    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updatePopulations(ModelViews& views,
                                                                        double const& dt,
                                                                        UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updatePopulations");
    resetMoments(ions);
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(ions, em, layout);
    else
        updateAndDepositAll_(ions, em, layout);
}




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositDomain_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositDomain_");

    Pusher_t pusher{dt_};

    auto& pool = PHARE::core::gpu::AsynchLauncherPool::INSTANCE();

    pool.check(views);
    pusher(views);

    // for (auto& pop : ions)
    // {
    //     auto& domain = pop.domainParticles();
    //     Pusher_t pusher; // dt? {layout, layout.meshSize(), dt_, pop.mass()};

    //     pusher.template push<ParticleType::Domain>(domain, em);
    //     interpolator_.particleToMesh(*domain, layout, pop.density(), pop.flux());

    //     auto pushGhosts = [&](auto& inputArray, bool copyInDomain = false) mutable {
    //         pusher.template push<ParticleType::Ghost>(tmp_particles_.replace_from(inputArray),
    //         em); interpolator_.particleToMesh(*tmp_particles_, layout, pop.density(),
    //         pop.flux()); if (copyInDomain)
    //             ParticleArrayService::copy_ghost_into_domain(tmp_particles_, domain);
    //     };

    //     pushGhosts(pop.patchGhostParticles(), true);
    //     pushGhosts(pop.levelGhostParticles());
    // }
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
/**
 * @brief IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositAll_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositAll_");


    Pusher_t pusher{dt_};

    pusher(views);

    // for (auto& pop : ions)
    // {
    //     auto& domain = pop.domainParticles();
    //     Pusher_t pusher{layout, layout.meshSize(), dt_, pop.mass()};
    //     pusher.template push<ParticleType::Domain>(domain, em);

    //     auto pushGhosts = [&](auto& particles) {
    //         pusher.template push<ParticleType::Ghost>(particles, em);
    //         ParticleArrayService::copy_ghost_into_domain(particles, domain);
    //     };
    //     pushGhosts(pop.patchGhostParticles());
    //     pushGhosts(pop.levelGhostParticles());

    //     interpolator_.particleToMesh(*domain, layout, pop.density(), pop.flux());
    // }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
