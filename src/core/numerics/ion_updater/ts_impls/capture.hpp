// IWYU pragma: private, include "core/numerics/ion_updater/ts_impls/entry.hpp"

#ifndef PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP
#define PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP


namespace PHARE::core::mkn
{


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdaterMultiTS
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
    IonUpdaterMultiTS(PHARE::initializer::PHAREDict const& /*dict*/) {}

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
        // tmp_particles_ = std::move(make_particles<Particles>(Box_t{}));
    }


private:
    template<typename ModelViews>
    void updateAndDepositDomain_(ModelViews& views);

    template<typename ModelViews>
    void updateAndDepositAll_(ModelViews& views);



    // dealloced on regridding/load balancing coarsest
    // Particles tmp_particles_;

    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updatePopulations(ModelViews& views,
                                                                        double const& dt,
                                                                        UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updatePopulations");

    for (auto& view : views)
        resetMoments(view.ions);
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(views);
    else
        updateAndDepositAll_(views);
}




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositDomain_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositDomain_");

    // constexpr static std::uint8_t N_ARRAYS       = 2;
    // constexpr static std::uint8_t DOMAIN_ID      = 0;
    // constexpr static std::uint8_t LEVEL_GHOST_ID = 2;

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU

    // std::uint16_t const group_size = N_ARRAYS * views[0].ions->size();
    // MultiBoris<ModelViews> in{dt_, views};
    // Pusher_t::move(in);

    // // finished moving particles on patch
    // in.streamer.group_barrier(group_size);

    // // add new domain particles
    // in.streamer.host_group_mutex(group_size, [&](auto const i) {
    //     auto is_domain_particles = in.particle_type[i] == DOMAIN_ID;
    //     if (is_domain_particles || in.particles[i]->size() == 0)
    //         return;

    //     auto copy_in = [&](auto const j) {
    //         ParticleArrayService::copy_ghost_into_domain(*in.particles[i], *in.particles[j]);
    //     };

    //     if (in.particle_type[i - 1] == DOMAIN_ID)
    //         copy_in(i - 1);
    //     else
    //         return;
    //     in.particles[i]->clear();
    //     in.pviews[i].clear();
    // });


    // auto pps     = in.pviews.data();
    // auto rhos    = in.rhos.data();
    // auto fluxes  = in.fluxes.data();
    // auto layouts = in.layouts.data();

    // // finished adding new domain particles
    // in.streamer.group_barrier(group_size);

    // in.streamer.async_dev_idx(N_ARRAYS, DOMAIN_ID, [=](auto const i) { // 0 = domain
    //     Interpolating_t::box_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
    // });
    // // no patch ghost as they're injected into domain
    // in.streamer.async_dev_idx(N_ARRAYS, LEVEL_GHOST_ID, [=](auto const i) { // 2 = level
    // ghosts
    //     Interpolating_t::box_kernel(pps[i], layouts[i], fluxes[i], rhos[i]);
    // });

    // in.streamer();
    // in.streamer.sync();
    // in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain_.txt");

#else
    // throw std::runtime_error("No available implementation")
#endif
    //
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositAll_(ModelViews& views)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositAll_");

    // constexpr static std::uint8_t N_ARRAYS  = 3;
    // constexpr static std::uint8_t DOMAIN_ID = 0;

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU

    using field_type = Ions::field_type;

    MultiBoris<ModelViews> in{dt_, views};
    Pusher_t::move(in);

    in.streamer.host([&](auto const i) mutable {
        auto& view = in.views[i];

        auto per_ghost = [&](auto& domain, auto& ghost) {
            ParticleArrayService::copy_ghost_into_domain(ghost, domain, view.layout);
            ghost.clear();
        };

        auto per_ghosts = [&](auto& domain, auto&... ghosts) { (per_ghost(domain, ghosts), ...); };

        for (auto& pop : view.ions)
            per_ghosts(pop.domainParticles(), pop.patchGhostParticles(), pop.levelGhostParticles());
    });


    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        in.streamer.host([&](auto const i) mutable {
            PHARE_LOG_LINE_SS("");
            auto& view     = in.views[i];
            using Launcher = gpu::ChunkLauncher<false>;

            for (auto& pop : view.ions)
            {
                auto const pps = *pop.domainParticles();
                PHARE_LOG_LINE_SS(pop.domainParticles().size() << " " << pps.size());
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = pps().size();
                launcher.ds  = 9 * 9 * 9 * 4 * 8; // !!!!
                launcher.stream(in.streamer.streams[i],
                                [=, layout = view.layout] __device__() mutable {
                                    Interpolating_t::template on_tiles<field_type>(pps, layout);
                                });
                in.streamer.streams[i].sync();
            }
        });
    }
    else // if CPU
    {
        in.streamer.host([&](auto const i) mutable {
            PHARE_LOG_LINE_SS("");
            auto& view = in.views[i];

            Interpolating_t interp;
            for (auto& pop : view.ions)
                interp.particleToMesh(pop.domainParticles(), view.layout, pop.density(),
                                      pop.flux());
        });
    }

    in.streamer().sync();
    // in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll_"
    //                        + std::string{Particles::type_id} + ".txt");

    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        for (auto& stream : in.streamer.streams)
            stream.sync();

        for (auto& view : views)
            for (auto& pop : view.ions)
                Interpolating_t::ts_reducer(pop.domainParticles(), view.layout, pop.flux(),
                                            pop.density());

        in.streamer().sync();
    }

#else
    // throw std::runtime_error("No available implementation")
#endif
    //
}


} // namespace PHARE::core::mkn


#endif // PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP
