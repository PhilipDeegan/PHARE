// IWYU pragma: private, include "core/numerics/ion_updater/ts_impls/entry.hpp"

#ifndef PHARE_ION_UPDATER_MULTI_TS_IMPLS_GPU_STATE_HPP
#define PHARE_ION_UPDATER_MULTI_TS_IMPLS_GPU_STATE_HPP


#include "core/def/phare_config.hpp"
#include "core/numerics/pusher/multi_boris.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/interpolator/interpolating.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"


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
        ions.computeChargeDensity();
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

    constexpr static std::uint8_t N_ARRAYS  = 3;
    constexpr static std::uint8_t DOMAIN_ID = 0;

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    using field_type = Ions::field_type;


    MultiBoris<ModelViews> in{dt_, views};
    Pusher_t::move(in);

    std::uint16_t const group_size = N_ARRAYS * views[0].ions.size();

    // finished moving particles on patch
    in.streamer.group_barrier(group_size);

    // add new domain particles
    in.streamer.host_group_mutex(group_size, [&](auto const i) {
        in.pviews[i].reset(*in.particles[i]);

        auto is_domain_particles = in.particle_type[i] == DOMAIN_ID;
        if (is_domain_particles || in.particles[i]->size() == 0)
            return;

        auto copy_in = [&](auto const j) {
            PHARE_LOG_LINE_SS(in.particles[j]->size());
            ParticleArrayService::copy_ghost_into_domain(*in.particles[i], *in.particles[j],
                                                         in.layouts[i]);
            in.pviews[j].reset(*in.particles[j]);
            PHARE_LOG_LINE_SS(in.particles[j]->size());
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

    // assert(GridLayout::nbrGhosts() == 2);
    // static_assert(GridLayout::nbrGhosts() == 2);

    std::uint8_t constexpr static INTERP_GPU_IMPL = 2;

    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED)
                  and INTERP_GPU_IMPL == 0)
    {
        auto const interper = [=] _PHARE_DEV_FN_(auto const i) mutable {
            Interpolating_t::chunk_kernel_ts(pps[i], layouts[i], fluxes[i], rhos[i]);
        };
        // int const maxbytes   = 65536; // 64 KB
        // auto static func_    = &mkn::gpu::global_d_kernel<decltype(interper), std::uint32_t
        // const>; void* const gpu_func = *reinterpret_cast<void**>(&func_);
        // hipFuncSetAttribute(gpu_func, hipFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
        in.streamer.async_dev_chunk_idx(N_ARRAYS, DOMAIN_ID, interper, 4, 9 * 9 * 9 * 1 * 8);
    }
    else if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED)
                       and INTERP_GPU_IMPL == 1)
    {
        // auto const interper = [=] _PHARE_DEV_FN_(auto const i) mutable {
        //     Interpolating_t::chunk_kernel_ts_all(pps[i], layouts[i], fluxes[i], rhos[i]);
        // };
        // in.streamer.async_dev_chunk_idx(N_ARRAYS, DOMAIN_ID, interper, 4, 9 * 9 * 9 * 4 * 8);
    }
    else if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED)
                       and INTERP_GPU_IMPL == 2)
    {
        in.streamer.async_dev_tiles(
            [=] _PHARE_DEV_FN_(auto const i) mutable {
                Interpolating_t::template on_tiles<field_type>(pps[i], layouts[i]);
            },
            9 * 9 * 9 * 4 * 8);
    }
    else
    {
        in.streamer.host([&](auto const i) {
            Interpolating_t interp;
            if (pps[i].size())
                interp.particleToMesh(pps[i], layouts[i], rhos[i], fluxes[i]);
        });
    }


    in.streamer().sync();

    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll_"
                           + std::string{Particles::type_id} + ".txt");



#else
    // throw std::runtime_error("No available implementation")
#endif
    //
}


} // namespace PHARE::core::mkn


#endif // PHARE_ION_UPDATER_MULTI_TS_IMPLS_GPU_STATE_HPP
