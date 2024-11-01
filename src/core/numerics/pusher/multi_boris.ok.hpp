#ifndef PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP
#define PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include <cstdint>
#include <cassert>
#include <stdio.h>
#include <iostream>

#include "core/errors.hpp"
#include "core/logger.hpp"

#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"

#include "core/numerics/pusher/boris_simpler.hpp"

#if PHARE_HAVE_MKN_GPU
#include "core/data/mkn.gpu.hpp"
#endif // PHARE_HAVE_MKN_GPU

#include "core/data/particles/arrays/particle_array_soa.hpp"

namespace PHARE::core
{

// !!naive!!
template<typename ModelViews>
struct MultiBoris
{
    using ModelView    = ModelViews::value_type;
    using GridLayout_t = ModelView::GridLayout_t;

#if PHARE_HAVE_MKN_GPU
    using ParticleArray_t     = ModelView::ParticleArray_t;
    static constexpr auto dim = ParticleArray_t::dimension;
    using Electromag_t        = ModelView::Electromag_t;
    using Vecfield_t          = Electromag_t::vecfield_type;
    using Field_t             = Vecfield_t::field_type;
    using ParticleArray_v     = typename ParticleArray_t::view_t;
    using Box_t               = Box<std::uint32_t, dim>;
    using Boxes_t             = std::vector<Box_t>;
    // = decltype(*(*views[0].ions).getRunTimeResourcesViewList()[0].domainParticles());
    using Particles_ptrs = std::vector<ParticleArray_t*>;
    // using StreamLauncher = gpu::BoxStreamLauncher<Boxes_t, Particles_ptrs>;
    using StreamLauncher = gpu::ThreadedBoxStreamLauncher<Boxes_t, Particles_ptrs>;
    // using StreamLauncher = mkn::gpu::StreamLauncher<>;

    static auto _particles(ModelViews& views)
    {
        std::vector<ParticleArray_t*> ptrs;
        auto add = [&](auto& ps) {
            // if (ps.size())
            ptrs.emplace_back(&ps);
        };
        auto all = [&](auto&... ps) { (add(ps), ...); };

        for (auto& view : views)
            for (auto& pop : *view.ions)
                all(pop.domainParticles(), pop.patchGhostParticles(), pop.levelGhostParticles());

        return ptrs;
    }

    MultiBoris(double const dt, ModelViews& _views)
        : views{_views}
        , particles{_particles(views)}
    {
        for (auto& view : views)
        {
            auto each = [&](auto const& type, auto& ps, auto& pop) mutable {
                // if (!ps.size())
                // return;
                particle_type.emplace_back(type);
                auto& domainview = pviews.emplace_back(*ps);
                boxes.emplace_back(domainview.local_box());
                layouts.emplace_back(view.layout);
                ems.emplace_back(*view.em);
                dto2ms.emplace_back(0.5 * dt / pop.mass());
                halfdt.emplace_back(mesh(view.layout.meshSize(), dt));
                rhos.emplace_back(pop.density());
                fluxes.emplace_back(pop.flux());
            };

            for (auto& pop : *view.ions)
            {
                each(0, pop.domainParticles(), pop);
                each(1, pop.patchGhostParticles(), pop);
                each(2, pop.levelGhostParticles(), pop);
            }
        }
    }

    void reset() {}

    double dt;
    ModelViews& views;

    // on host
    Particles_ptrs particles;
    std::vector<std::uint16_t> particle_type;
    Boxes_t boxes;

    // on device
    gpu::Vec_t<ParticleArray_v> pviews;
    gpu::Vec_t<GridLayout_t> layouts;
    gpu::Vec_t<Electromag_t> ems;
    gpu::Vec_t<Field_t> rhos;
    gpu::Vec_t<Vecfield_t> fluxes;
    gpu::Vec_t<std::array<double, dim>> halfdt;
    gpu::Vec_t<double> dto2ms;

    StreamLauncher streamer{particles, boxes, 5};

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }

#endif // PHARE_HAVE_MKN_GPU
};

template<typename GridLayout, typename Particles, typename Electromag, typename Interpolator>
class MultiBorisPusher
{
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;
    using This                = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator>;

public:
    template<typename ModelViews>
    static void move(MultiBoris<ModelViews>& in)
    {
#if PHARE_HAVE_MKN_GPU

        using ParticleArray_v            = typename MultiBoris<ModelViews>::ParticleArray_v;
        using Box_t                      = Box<std::uint32_t, dim>;
        using StreamLauncher             = typename MultiBoris<ModelViews>::StreamLauncher;
        using BoxLauncher                = gpu::BoxCellNLauncher<Box_t, false>;
        static constexpr auto alloc_mode = Particles_t::alloc_mode;

        auto const emps         = in.ems.data();
        auto const dto2mspp     = in.dto2ms.data();
        auto const layoutps     = in.layouts.data();
        auto const halfDtOverDl = in.halfdt.data();
        auto pps                = in.pviews.data();


        auto per_particle = [=] _PHARE_DEV_FN_(auto const& i) mutable {
            auto const blockidx  = BoxLauncher::block_idx();
            auto const threadIdx = BoxLauncher::thread_idx();
            auto const& dto2m    = dto2mspp[i];
            auto const& layout   = layoutps[i];
            auto& view           = pps[i];
            auto const& lobox    = view.local_box();
            auto const& locell   = *(lobox.begin() + blockidx);
            auto& parts          = view(locell);
            if (threadIdx >= parts.size())
                return;

            if constexpr (any_in(ParticleArray_v::layout_mode, LayoutMode::SoA, LayoutMode::SoAPC))
            {
                detail::SoAZipParticle particle{parts, threadIdx};
                auto const og_iCell = particle.iCell();

                particle.iCell() = advancePosition_<alloc_mode>(particle, halfDtOverDl[i]);
                {
                    Interpolator interp;
                    boris_accelerate(particle, interp.m2p(particle, emps[i], layout), dto2m);
                }

                particle.iCell() = advancePosition_<alloc_mode>(particle, halfDtOverDl[i]);

                if (!array_equals(particle.iCell(), og_iCell))
                    view.icell_changer(particle, locell, threadIdx, particle.iCell());
            }
            else
            {
                auto& particle      = parts[threadIdx];
                auto const og_iCell = particle.iCell();

                particle.iCell() = advancePosition_<alloc_mode>(particle, halfDtOverDl[i]);
                {
                    Interpolator interp;
                    boris_accelerate(particle, interp.m2p(particle, emps[i], layout), dto2m);
                }

                particle.iCell() = advancePosition_<alloc_mode>(particle, halfDtOverDl[i]);

                if (!array_equals(particle.iCell(), og_iCell))
                    view.icell_changer(particle, locell, threadIdx, particle.iCell());
            }
        };


        auto& streamer = in.streamer;
        auto ip        = &in; // used in lambdas, copy address! NO REF!

        streamer.async_dev([=] _PHARE_DEV_FN_(auto const i) mutable { per_particle(i); });
        streamer.host([ip = ip](auto const i) mutable {
            constexpr static std::uint32_t PHASE = 1;

            if (ip->particles[i]->size() == 0)
                return;

            if (ip->particle_type[i] == 0) //  :(
                This::template sync<PHASE, ParticleType::Domain>(ip, i);
            else if (any_in(ip->particle_type[i], 1, 2))
                This::template sync<PHASE, ParticleType::Ghost>(ip, i);
            else
                throw std::runtime_error("No impl");
        });

#endif
    }

    template<auto phase, auto type>
    static void sync(auto p, auto const& i)
    {
#if PHARE_HAVE_MKN_GPU
        assert(i < p->pviews.size() and i < p->streamer.streams.size());
        p->pviews[i].template sync<phase, type>(p->streamer.streams[i]);
        ParticleArrayService::sync<phase, type>(*p->particles[i]);
#endif
    }
};

} // namespace PHARE::core


#endif /* PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP */
