#ifndef PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP
#define PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"

#include "core/numerics/pusher/boris/move.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_service.hpp"


#include <cmath>
#include <cassert>
#include <cstddef>


namespace PHARE::core
{


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
        using ParticleArray_v            = typename MultiBoris<ModelViews>::ParticleArray_v;
        static constexpr auto alloc_mode = Particles_t::alloc_mode;

        auto const emps         = in.ems.data();
        auto const dto2mspp     = in.dto2ms.data();
        auto const layoutps     = in.layouts.data();
        auto const halfDtOverDl = in.halfdt.data();
        auto pps                = in.pviews.data();

        auto per_particle = [=] _PHARE_ALL_FN_(auto&&... args) mutable {
            auto const& [particle, view, tile_cell, pidx, i] = std::forward_as_tuple(args...);

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
            {
                Interpolator interp;
                boris::accelerate(particle, interp.m2p(particle, emps[i], layoutps[i]),
                                  dto2mspp[i]);
            }
            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

            auto const new_cell   = view.local_tile_cell(particle.iCell());
            bool const moved_tile = !array_equals(new_cell, tile_cell);
            if (moved_tile)
                view.icell_changer(particle, tile_cell, pidx, particle.iCell());
        };

        auto per_any_particle = [=] _PHARE_DEV_FN_(auto& particles, auto&&... args) mutable {
            auto const& pidx = std::get<2>(std::forward_as_tuple(args...));

            using enum LayoutMode;
            if constexpr (any_in(ParticleArray_v::layout_mode, SoA, SoAPC, SoATS))
                per_particle(detail::SoAZipParticle{particles, pidx}, args...);
            else
                per_particle(particles[pidx], args...);
        };

        auto per_tile_chunk = [=] _PHARE_DEV_FN_(auto const& i) mutable {
            Point<std::uint32_t, 3> const rcell{threadIdx.z, threadIdx.y, threadIdx.x};
            auto& view        = pps[i];
            auto& tile        = *view().at(rcell + 2);
            auto& parts       = tile();
            auto const tile_p = rcell % 2;
            auto const tidx   = tile_p[2] + tile_p[1] * 2 + tile_p[0] * 2 * 2;
            auto const each   = parts.size() / 8;
            for (std::size_t pid = 0; pid < each; ++pid)
                per_any_particle(parts, view, view.local_cell(tile.lower), pid * 8 + tidx, i);
        };

        auto per_tile = [=] _PHARE_DEV_FN_(auto const& i) mutable {
            auto& view      = pps[i];
            auto& tile      = view()[blockIdx.x];
            auto& parts     = tile();
            auto const tidx = threadIdx.x;
            auto const ws   = kernel::warp_size();
            auto const each = parts.size() / ws;

            std::size_t pid = 0;
            for (; pid < each; ++pid)
                per_any_particle(parts, view, view.local_cell(tile.lower), pid * ws + tidx, i);
            if (tidx < parts.size() - (ws * each))
                per_any_particle(parts, view, view.local_cell(tile.lower), pid * ws + tidx, i);
        };


        auto& streamer = in.streamer;
        auto ip        = &in; // used in lambdas, copy address! NO REF!

        if constexpr (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            // streamer.async_dev_chunk([=] _PHARE_DEV_FN_(auto const i) mutable { per_tile(i); },
            // 4);
            streamer.async_dev_tiles([=] _PHARE_DEV_FN_(auto const i) mutable { per_tile(i); });

            streamer.host([ip = ip](auto const i) mutable {
                if (ip->particles[i]->size() == 0)
                    return;

                constexpr static std::uint32_t PHASE = 1;
                if (ip->particle_type[i] == 0) //  :(
                    This::template sync<PHASE, ParticleType::Domain>(ip, i);
                else if (any_in(ip->particle_type[i], 1, 2))
                    This::template sync<PHASE, ParticleType::Ghost>(ip, i);
                else
                    throw std::runtime_error("No impl");
            });
        }
        else
        {
            // do cpu
            streamer.host([=](auto const i) mutable {
                if (ip->particles[i]->size() == 0)
                    return;

                auto& particles_ts = *ip->particles[i];

                for (auto& tile : particles_ts())
                {
                    auto const locell  = particles_ts.local_cell(tile.lower);
                    auto const& nparts = particles_ts.size(locell);
                    for (std::size_t pidx = 0; pidx < nparts; ++pidx)
                        per_any_particle(tile(), particles_ts, locell, pidx, i);
                }
            });

            // separation for profiling the same as gpu
            streamer.host([=](auto const i) mutable {
                constexpr static std::uint32_t PHASE = 1;

                auto& particles_ts = *ip->particles[i];
                if (ip->particle_type[i] == 0)
                    ParticleArrayService::sync<PHASE, ParticleType::Domain>(particles_ts);
                else if (any_in(ip->particle_type[i], 1, 2))
                    ParticleArrayService::sync<PHASE, ParticleType::Ghost>(particles_ts);
                else
                    throw std::runtime_error("No impl");
            });
        }
    }

    template<auto phase, auto type>
    static void sync(auto p, auto const& i)
    {
        assert(i < p->pviews.size() and i < p->streamer.streams.size());
        p->pviews[i].template sync<phase, type>(p->streamer.streams[i]);
        ParticleArrayService::sync<phase, type>(*p->particles[i]);
    }
};

} // namespace PHARE::core


#endif /* PHARE_CORE_PUSHER_MULTI_BORIS_2_HPP */
