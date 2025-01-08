#ifndef PHARE_ION_UPDATER_MULTI_PC_HPP
#define PHARE_ION_UPDATER_MULTI_PC_HPP

#include "ion_updater_def.hpp"

#include "core/numerics/pusher/multi_boris.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

#if PHARE_HAVE_MKN_GPU || defined(_MKN_WITH_MKN_KUL_)
#include "mkn/kul/os.hpp"
#endif // PHARE_HAVE_MKN_GPU


namespace PHARE::core::detail
{
auto static const timings_dir_str
    = get_env_as("PHARE_ASYNC_TIMES", std::string{".phare/async/multi_updater"});

static bool ion_updater_pc_setup = []() {
    mkn::kul::Dir timings{timings_dir_str};
    timings.mk();
    return true;
}();

} // namespace PHARE::core::detail


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

    constexpr static std::uint8_t N_ARRAYS       = 1;
    constexpr static std::uint8_t DOMAIN_ID      = 0;
    constexpr static std::uint8_t LEVEL_GHOST_ID = 2;

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    std::uint16_t const group_size = N_ARRAYS * views[0].ions->size();
    MultiBoris<ModelViews> in{dt_, views};
    Pusher_t::move(in);

    // finished moving particles on patch
    in.streamer.group_barrier(group_size);

    // add new domain particles
    in.streamer.host_group_mutex(group_size, [&](auto const i) {
        auto is_domain_particles = in.particle_type[i] == DOMAIN_ID;
        if (is_domain_particles || in.particles[i]->size() == 0)
            return;

        auto copy_in = [&](auto const j) {
            ParticleArrayService::copy_ghost_into_domain(*in.particles[i], *in.particles[j]);
        };

        if (in.particle_type[i - 1] == DOMAIN_ID)
            copy_in(i - 1);
        else
            return;
        in.particles[i]->clear();
        in.pviews[i].clear();
    });


    auto pps     = in.pviews.data();
    auto rhos    = in.rhos.data();
    auto fluxes  = in.fluxes.data();
    auto layouts = in.layouts.data();

    // finished adding new domain particles
    in.streamer.group_barrier(group_size);

    in.streamer.async_dev_idx(N_ARRAYS, DOMAIN_ID, [=](auto const i) { // 0 = domain
        Interpolating_t::box_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
    });
    // no patch ghost as they're injected into domain
    in.streamer.async_dev_idx(N_ARRAYS, LEVEL_GHOST_ID, [=](auto const i) { // 2 = level ghosts
        Interpolating_t::box_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
    });

    in.streamer();
    in.streamer.sync();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain_.txt");

#else
        // throw std::runtime_error("No available implementation")
#endif
    //
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositAll_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositAll_");

    constexpr static std::uint8_t N_ARRAYS  = 1;
    constexpr static std::uint8_t DOMAIN_ID = 0;
    // constexpr static std::uint8_t LEVEL_GHOST_ID = 2;

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    MultiBoris<ModelViews> in{dt_, views};
    Pusher_t::move(in);

    std::uint16_t const group_size = N_ARRAYS * views[0].ions->size();

    // finished moving particles on patch
    in.streamer.group_barrier(group_size);

    // add new domain particles
    in.streamer.host_group_mutex(group_size, [&](auto const i) {
        // hipStreamSynchronize(in.streamer.streams[i]);

        auto is_domain_particles = in.particle_type[i] == DOMAIN_ID;
        if (is_domain_particles || in.particles[i]->size() == 0)
            return;

        auto copy_in = [&](auto const j) {
            ParticleArrayService::copy_ghost_into_domain(*in.particles[i], *in.particles[j]);
            in.pviews[j].reset(*in.particles[j]);
        };

        if (in.particle_type[i - 1] == DOMAIN_ID)
            copy_in(i - 1);
        else if (in.particle_type[i - 2] == DOMAIN_ID)
            copy_in(i - 2);

        // means there is no particle to mesh for ghosts
        in.particles[i]->clear();
        in.pviews[i].clear();
    });


    auto pps     = in.pviews.data();
    auto rhos    = in.rhos.data();
    auto fluxes  = in.fluxes.data();
    auto layouts = in.layouts.data();

    // finished adding new domain particles
    in.streamer.group_barrier(group_size);

    assert(GridLayout::nbrGhosts() == 2);
    static_assert(GridLayout::nbrGhosts() == 2);

    in.streamer.async_dev_chunk(
        [=] _PHARE_DEV_FN_(auto const i) mutable {
            Interpolating_t::chunk_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
        },
        4, 9 * 9 * 9 * 4 * 8);

    in.streamer();
    in.streamer.sync();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll_"
                           + std::string{Particles::type_id} + ".txt");

#else
        // throw std::runtime_error("No available implementation")
#endif
    //
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
