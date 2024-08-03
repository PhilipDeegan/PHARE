#ifndef PHARE_ION_UPDATER_MULTI_PC_HPP
#define PHARE_ION_UPDATER_MULTI_PC_HPP

#include "ion_updater_def.hpp"

#include "core/numerics/pusher/multi_boris.hpp"
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

    using Pusher_t = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator_t>;

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
    // resetMoments(ions);
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(views);
    else
        updateAndDepositAll_(views);
}




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositDomain_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositDomain_");

    // Pusher_t pusher{dt_};
    // auto launcher = pusher.move_async(views);

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
    static_assert(PHARE_HAVE_MKN_GPU);
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositAll_");

#if PHARE_HAVE_MKN_GPU
    MultiBoris<ModelViews> in{dt_, views};
    Pusher_t::move(in);

    std::mutex mute;

    in.streamer.host([&](auto const i) {
        auto is_domain_particles = in.particle_type[i] == 0;
        if (is_domain_particles || in.particles[i]->size() == 0)
            return;

        {
            std::scoped_lock<std::mutex> lk(mute); // copy back to domain

            auto copy_in = [&](auto const j) {
                ParticleArrayService::copy_ghost_into_domain(*in.particles[i], *in.particles[j]);
            };

            // do better
            if (in.particle_type[i - 1] == 0)
                copy_in(i - 1);
            else if (in.particle_type[i - 2] == 0)
                copy_in(i - 2);
        }

        in.particles[i]->clear();
        in.pviews[i].clear();
    });


    auto pps     = in.pviews.data();
    auto rhos    = in.rhos.data();
    auto fluxes  = in.fluxes.data();
    auto layouts = in.layouts.data();

    in.streamer.barrier(); // barrier_group for one patch?

    in.streamer.async_dev([=](auto const i) { //
        Interpolating_t::box_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
    });

    in.streamer();
    in.streamer.sync();

#else
    throw std::runtime_error("No available implementation")
#endif
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
