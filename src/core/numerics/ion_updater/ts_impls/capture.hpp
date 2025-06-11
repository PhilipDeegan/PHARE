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
    // static_assert(Particles::alloc_mode == AllocatorMode::GPU_UNIFIED);

public:
    auto static constexpr dimension    = GridLayout::dimension;
    auto static constexpr interp_order = GridLayout::interp_order;
    using Box_t                        = Box<int, dimension>;
    using ParticleArray_t              = Particles;
    using Interpolator_t               = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t              = Interpolating<Particles, interp_order, atomic_interp>;

    using Pusher_t = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator_t>;

private:
    Interpolating_t interpolator_;


public:
    IonUpdaterMultiTS(PHARE::initializer::PHAREDict const& /*dict*/) {}

    template<typename ModelViews, typename Boxing_t>
    void updatePopulations(ModelViews&, Boxing_t const&, double const&,
                           UpdaterMode = UpdaterMode::all);


    static void updateIons(Ions& ions)
    {
        // fixMomentGhosts(ions, layout);
        ions.computeDensity();
        ions.computeBulkVelocity();
    }

    void reset() {} // noop


private:
    template<typename ModelViews, typename Boxing_t>
    void updateAndDepositDomain_(ModelViews& views, Boxing_t const&);

    template<typename ModelViews, typename Boxing_t>
    void updateAndDepositAll_(ModelViews& views, Boxing_t const&);


    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updatePopulations(ModelViews& views,
                                                                        Boxing_t const& boxings,
                                                                        double const& dt,
                                                                        UpdaterMode mode)
{
    PHARE_LOG_SCOPE(2, "IonUpdaterMultiTS::updatePopulations");

    for (auto& view : views)
        resetMoments(view.ions);
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(views, boxings);
    else
        updateAndDepositAll_(views, boxings);
}




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositDomain_(
    ModelViews& views, Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositDomain_");

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    using MultiBoris_t = MultiBoris<ModelViews>;
    MultiBoris<ModelViews> in{dt_, views};

    Pusher_t::template move<MultiBorisMode::COPY>(in);

    in.streamer.host([&](auto const i) mutable {
        auto& view           = in.views[i];
        auto const& patch_id = view.patchID();
        assert(boxings.count(patch_id));
        auto const& patch_boxings = boxings.at(patch_id);

        for (std::size_t j = 0; j < view.ions.size(); ++j)
        {
            auto const popidx = view.ions.size() * i + j;
            delete_particles<false>(MultiBoris_t::domains[popidx], patch_boxings.nonLevelGhostBox);

            move_in_domain(MultiBoris_t::domains[popidx], MultiBoris_t::levelGhosts[popidx],
                           patch_boxings.domainBox);
        }
    });

    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        throw std::runtime_error("finish");
    }
    else
    {
        in.streamer.host([&](auto const i) mutable {
            PHARE_LOG_LINE_SS("");
            auto& view = in.views[i];

            std::size_t interped = 0;
            Interpolating_t interp;

            for (std::size_t j = 0; j < view.ions.size(); ++j)
            {
                auto const popidx = view.ions.size() * i + j;
                auto& pop         = view.ions[j];
                interp.particleToMesh(MultiBoris_t::domains[popidx], view.layout, pop.density(),
                                      pop.flux());

                interped += MultiBoris_t::domains[popidx].size();
            }
            PHARE_LOG_LINE_SS(interped);
        });
    }

    in.streamer.sync();
    // in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain_.txt");

#else
    throw std::runtime_error("No available implementation")
#endif
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelViews, typename Boxing_t>
void IonUpdaterMultiTS<Ions, Electromag, GridLayout>::updateAndDepositAll_(ModelViews& views,
                                                                           Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiTS::updateAndDepositAll_");

    if (views.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU

    using field_type = Ions::field_type::value_type;


    MultiBoris<ModelViews> in{dt_, views};

    Pusher_t::move(in);

    in.streamer.host([&](auto const i) mutable {
        auto& view           = in.views[i];
        auto const& patch_id = view.patchID();
        assert(boxings.count(patch_id));
        auto const& patch_boxings = boxings.at(patch_id);

        auto const per_pop = [&](auto& pop) {
            auto& domain = pop.domainParticles();

            assert(domain.size());
            delete_particles<false>(domain, patch_boxings.nonLevelGhostBox);

            move_in_ghost_layer(pop.patchGhostParticles(), domain, patch_boxings.domainBox,
                                patch_boxings.nonLevelGhostBox);

            move_in_domain(domain, pop.levelGhostParticles(), patch_boxings.domainBox);

            delete_particles<false>(pop.levelGhostParticles(), patch_boxings.ghostBox);
            delete_particles<false>(domain, patch_boxings.domainBox);

            assert(domain.size());

            for (auto const& tile : domain())
                for (auto const& p : tile())
                {
                    assert(isIn(p, patch_boxings.domainBox));
                    assert(isIn(p, tile));
                }
            for (auto const& tile : pop.levelGhostParticles()())
                for (auto const& p : tile())
                {
                    assert(isIn(p, patch_boxings.ghostBox)
                           and not isIn(p, patch_boxings.domainBox));
                    assert(not isIn(p, tile));
                }
            for (auto const& tile : pop.patchGhostParticles()())
                for (auto const& p : tile())
                {
                    assert(isIn(p, patch_boxings.ghostBox)
                           and not isIn(p, patch_boxings.domainBox));
                    assert(not isIn(p, tile));
                }
        };

        for (auto& pop : view.ions)
            per_pop(pop);
    });


    if constexpr (any_in(Particles::alloc_mode, AllocatorMode::GPU_UNIFIED))
    {
        in.streamer.host([&](auto const i) mutable {
            PHARE_LOG_LINE_SS("");
            auto& view     = in.views[i];
            using Launcher = gpu::ChunkLauncher<false>;

            for (auto& pop : view.ions)
            {
                auto const domain = *pop.domainParticles();
                auto const pghost = *pop.patchGhostParticles();
                PHARE_LOG_LINE_SS(pop.domainParticles().size() << " " << domain.size());
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = domain().size();
                launcher.ds  = 9 * 9 * 9 * 4 * 8; // !!!! NO LONGER VALID non-homogenuous tile sizes
                launcher.stream(in.streamer.streams[i],
                                [=, layout = view.layout] __device__() mutable {
                                    Interpolating_t::template on_tiles<field_type>(domain, layout);
                                    Interpolating_t::template on_tiles<field_type>(pghost, layout);
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

            std::size_t interped = 0;
            Interpolating_t interp;
            for (auto& pop : view.ions)
            {
                interp.particleToMesh(pop.domainParticles(), view.layout, pop.density(),
                                      pop.flux());
                interp.particleToMesh(pop.patchGhostParticles(), view.layout, pop.density(),
                                      pop.flux());

                PHARE_LOG_LINE_SS(pop.domainParticles().size());
                PHARE_LOG_LINE_SS(pop.patchGhostParticles().size());

                interped += pop.domainParticles().size();
                interped += pop.patchGhostParticles().size();
            }

            PHARE_LOG_LINE_SS(interped);
        });
    }

    in.streamer.join(false);

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

    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain_.txt");

#else
    throw std::runtime_error("No available implementation")
#endif
    //
}


} // namespace PHARE::core::mkn


#endif // PHARE_ION_UPDATER_MULTI_TS_IMPLS_CAPTURE_HPP
