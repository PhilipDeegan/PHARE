// IWYU pragma: private, include "core/numerics/ion_updater/ts_impls/entry.hpp"

#ifndef PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP
#define PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP


#include "core/def/phare_config.hpp"
#include "core/numerics/pusher/multi_boris.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/interpolator/interpolating.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"
#include "core/data/particles/particle_array_exporter.hpp"


namespace PHARE::core::mkn
{


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdaterMultiTS
{
    using Particles                     = Ions::particle_array_type;
    bool constexpr static atomic_interp = Particles::alloc_mode == AllocatorMode::GPU_UNIFIED;

public:
    auto static constexpr dimension    = GridLayout::dimension;
    auto static constexpr interp_order = GridLayout::interp_order;
    using Box_t                        = Box<int, dimension>;
    using ParticleArray_t              = Particles;
    using Interpolator_t               = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t              = Interpolating<Particles, interp_order, atomic_interp>;
    using Pusher_t   = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator_t>;
    using Vecfield_t = Electromag::vecfield_type;
    using Field_t    = Vecfield_t::field_type;
    using Tile_vt    = Field_t::value_type;

private:
    Interpolating_t interpolator_;


public:
    IonUpdaterMultiTS(PHARE::initializer::PHAREDict const& /*dict*/) {}

    template<typename ModelAccessor_t, typename Boxing_t>
    void updatePopulations(ModelAccessor_t&, Boxing_t const&, double const&,
                           UpdaterMode = UpdaterMode::all);


    static void updateIons(Ions& ions)
    {
        ions.computeChargeDensity();
        ions.computeBulkVelocity();
    }

    void reset() {}


private:
    template<typename ModelAccessor_t, typename Boxing_t>
    void updateAndDepositDomain_(ModelAccessor_t&, Boxing_t const&);

    template<typename ModelAccessor_t, typename Boxing_t>
    void updateAndDepositAll_(ModelAccessor_t&, Boxing_t const&);


    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updatePopulations(ModelAccessor_t& accessor,
                                                                        Boxing_t const& boxings,
                                                                        double const& dt,
                                                                        UpdaterMode mode)
{
    PHARE_LOG_SCOPE(2, "IonUpdaterMultiTS::updatePopulations");

    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        auto view      = accessor[i];
        auto [ions, _] = view.args;
        resetMoments(ions);
    }
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(accessor, boxings);
    else
        updateAndDepositAll_(accessor, boxings);
}




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositDomain_(
    ModelAccessor_t& accessor, Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositDomain_");

    if (accessor.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    using MultiBoris_t = MultiBoris<ModelAccessor_t>;
    MultiBoris_t in{dt_, accessor};

    Pusher_t::template move<MultiBorisMode::COPY>(in, boxings);

    in.streamer.host([&](auto const i) mutable {
        if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
        {
            auto view                 = accessor[i];
            auto const patch_id       = view.patchID();
            auto const& patch_boxings = boxings.at(patch_id);
            auto [ions, _]            = view.args;

            for (std::size_t j = 0; j < ions.size(); ++j)
            {
                auto const popidx = ions.size() * i + j;
                delete_particles_not_in(*MultiBoris_t::domains[popidx],
                                        patch_boxings.nonLevelGhostBox);
                move_in_domain(*MultiBoris_t::domains[popidx], *MultiBoris_t::levelGhosts[popidx],
                               patch_boxings.domainBox);
            }
        }
    });

    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        in.streamer.host([&](auto const i) mutable {
            auto view      = accessor[i];
            auto [ions, _] = view.args;
            for (std::size_t j = 0; j < ions.size(); ++j)
            {
                auto const popidx      = ions.size() * i + j;
                auto& pop              = ions[j];
                std::uint32_t const ds = ions.chargeDensity().max_tile_size();
                auto pps               = **MultiBoris_t::domains[popidx];
                using Launcher         = gpu::ChunkLauncher<false>;
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = pps().size();
                launcher.ds  = ds * 5 * 8;
                assert(launcher.ds < 65000);
                auto rhop = pop.particleDensity();
                auto rhoc = pop.chargeDensity();
                auto flux = *pop.flux();
                launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                    Interpolating_t::template on_tiles<Tile_vt>(pps, flux, rhop, rhoc);
                });
            }
        });
    }
    else
    {
        in.streamer.host([&]([[maybe_unused]] auto const i) mutable {});
    }

    in.streamer.join();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain.txt");

#else
    throw std::runtime_error("No available implementation");
#endif
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositAll_(
    ModelAccessor_t& accessor, Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositAll_");

    if (accessor.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU

    MultiBoris<ModelAccessor_t> in{dt_, accessor};

    Pusher_t::move(in, boxings);

    in.streamer.host([&](auto const i) mutable {
        auto view                 = accessor[i];
        auto [ions, _]            = view.args;
        auto const patch_id       = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);

        auto const per_pop = [&](auto& pop) {
            auto& domain = pop.domainParticles();
            delete_particles_not_in(domain, patch_boxings.nonLevelGhostBox);
            move_in_ghost_layer(pop.patchGhostParticles(), domain, patch_boxings.domainBox,
                                patch_boxings.nonLevelGhostBox);
            move_in_domain(domain, pop.levelGhostParticles(), patch_boxings.domainBox);
            delete_particles_not_in(pop.levelGhostParticles(), patch_boxings.ghostBox);
            delete_particles_not_in(domain, patch_boxings.domainBox);
        };

        for (auto& pop : ions)
            per_pop(pop);
    });

    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        in.streamer.host([&](auto const i) mutable {
            auto view      = accessor[i];
            auto [ions, _] = view.args;
            for (std::size_t j = 0; j < ions.size(); ++j)
            {
                auto& pop              = ions[j];
                std::uint32_t const ds = ions.chargeDensity().max_tile_size();
                auto const domain      = *pop.domainParticles();
                auto const pghost      = *pop.patchGhostParticles();
                using Launcher         = gpu::ChunkLauncher<false>;
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = domain().size();
                launcher.ds  = ds * 5 * 8;
                assert(launcher.ds < 65000);
                auto pdensity = pop.particleDensity();
                auto cdensity = pop.chargeDensity();
                auto flux     = *pop.flux();
                launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                    Interpolating_t::template on_tiles<Tile_vt>(domain, flux, pdensity, cdensity);
                });
            }
        });
    }
    else
    {
        in.streamer.host([&](auto const i) mutable {
            auto view            = accessor[i];
            auto [ions, _]       = view.args;
            auto const patch_id  = view.patchID();
            auto const& boxing_i = boxings.at(patch_id);
            Interpolating_t interp;
            for (auto& pop : ions)
            {
                interp.particleToMesh(pop.domainParticles(), boxing_i.layout, pop.particleDensity(),
                                      pop.chargeDensity(), pop.flux());
                interp.particleToMesh(pop.patchGhostParticles(), boxing_i.layout,
                                      pop.particleDensity(), pop.chargeDensity(), pop.flux());
            }
        });
    }

    in.streamer.join();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll.txt");

#else
    throw std::runtime_error("No available implementation");
#endif
}


} // namespace PHARE::core::mkn


#endif // PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP
