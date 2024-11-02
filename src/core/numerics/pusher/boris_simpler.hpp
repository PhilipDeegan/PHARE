#ifndef PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP
#define PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP

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
#include "core/vector.hpp"


#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"

#include "core/numerics/pusher/pusher.hpp"


#if PHARE_HAVE_MKN_GPU
#include "core/data/mkn.gpu.hpp"
#endif // PHARE_HAVE_MKN_GPU

namespace PHARE::core
{


template<auto alloc, typename Particle_t, typename Float, std::size_t dim>
auto static advancePosition_(Particle_t& p,
                             std::array<Float, dim> const halfDtOverDl) _PHARE_ALL_FN_
{
    auto newCell = p.iCell();

    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        Float delta = p.delta()[iDim] + static_cast<Float>(halfDtOverDl[iDim] * p.v()[iDim]);
        if constexpr (alloc == AllocatorMode::CPU)
            if (std::abs(delta) > 2)
                throw_runtime_error("Error, particle moves more than 1 cell, delta >2");

        Float iCell     = std::floor(delta);
        p.delta()[iDim] = delta - iCell;
        newCell[iDim]   = static_cast<int>(iCell + p.iCell()[iDim]);
    }

    return newCell;
}


template<typename Particle, typename Float, typename ParticleEB>
void boris_accelerate(Particle& p, ParticleEB const& eb, Float const& dto2m) _PHARE_ALL_FN_
{
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    Float const coef1 = p.charge() * dto2m;

    // We now apply the 3 steps of the BORIS PUSHER
    // 1st half push of the electric field
    Float const velx1 = p.v()[0] + coef1 * pEx;
    Float const vely1 = p.v()[1] + coef1 * pEy;
    Float const velz1 = p.v()[2] + coef1 * pEz;

    // preparing variables for magnetic rotation
    Float const rx = coef1 * pBx;
    Float const ry = coef1 * pBy;
    Float const rz = coef1 * pBz;

    Float const rx2  = rx * rx;
    Float const ry2  = ry * ry;
    Float const rz2  = rz * rz;
    Float const rxry = rx * ry;
    Float const rxrz = rx * rz;
    Float const ryrz = ry * rz;

    Float const invDet = one / (one + rx2 + ry2 + rz2);

    // preparing rotation matrix due to the magnetic field
    // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
    Float const mxx = one + rx2 - ry2 - rz2;
    Float const mxy = two * (rxry + rz);
    Float const mxz = two * (rxrz - ry);

    Float const myx = two * (rxry - rz);
    Float const myy = one + ry2 - rx2 - rz2;
    Float const myz = two * (ryrz + rx);

    Float const mzx = two * (rxrz + ry);
    Float const mzy = two * (ryrz - rx);
    Float const mzz = one + rz2 - rx2 - ry2;

    // magnetic rotation
    Float const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
    Float const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
    Float const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;

    // 2nd half push of the electric field / Update particle velocity
    p.v()[0] = velx2 + coef1 * pEx;
    p.v()[1] = vely2 + coef1 * pEy;
    p.v()[2] = velz2 + coef1 * pEz;
}




template<std::size_t dim, typename Interpolator, typename GridLayout>
class SimpleBorisPusher
{
    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }

public:
    SimpleBorisPusher(                     //
        GridLayout const& layout,          //
        std::array<double, dim> const& ms, //
        double const& ts,                  //
        double const& mass)
        : layout_{layout}
        , halfDtOverDl_{mesh(ms, ts)}
        , mass_{mass}
        , dt_{ts}
    {
        PHARE_ASSERT(mass_ > 0);
        PHARE_ASSERT(dt_ > 0);
        PHARE_ASSERT(dto2m > 0);
    }


    template<typename Particles, typename Electromag>
    void operator()(Particles& particles, Electromag const& em) const
    {
        move<ParticleType::Domain>(particles, em);
    }

    template<auto type, typename Particles, typename Electromag>
    void push(Particles& particles, Electromag const& em) const
    {
        move<type>(particles, em);
    }



    template<auto type, typename Particles, typename Electromag>
    void move(Particles& particles, Electromag const& em) const _PHARE_ALL_FN_
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;
        // PHARE_LOG_LINE_STR(particles.size());
        if (particles.size() == 0)
            return; // noop

        PHARE_LOG_SCOPE(1, "boris::move");

        if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        {
            if constexpr (alloc_mode == AllocatorMode::CPU)
                move_aos_pc_cpu<type>(particles, em);
            else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
                move_aos_pc_gpu<type>(particles, em);
            else
                throw std::runtime_error("no");
        }
        else if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            for (auto& p : particles)
                pp<alloc_mode>(p, em, interpolator_, layout_, halfDtOverDl_, dto2m);
        }
        else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(            //
                auto view = *particles;    //
                PHARE_ASSERT(view.size()); //
                mkn::gpu::GDLauncher{particles.size()}([=, layout = layout_,
                                                        halfDtOverDl = halfDtOverDl_,
                                                        dto2m_ = dto2m] _PHARE_ALL_FN_() mutable {
                    Interpolator interp;
                    auto particle = view.begin() + mkn::gpu::idx();
                    pp<alloc_mode>(particle, em, interp, layout, halfDtOverDl, dto2m_);
                }); //
            )
        }
        else
            throw std::runtime_error("no");
    }


    template<auto alloc, typename P, typename E, typename I, typename G>
    void static pp(P& particle, E const& em, I& interp, G const& layout,
                   std::array<double, dim> const halfDtOverDl, double const& dto2m) _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_STR(particle.copy());
        particle.iCell() = advancePosition_<alloc>(particle, halfDtOverDl);
        boris_accelerate(particle, interp.m2p(particle, em, layout), dto2m);
        particle.iCell() = advancePosition_<alloc>(particle, halfDtOverDl);
        // PHARE_LOG_LINE_STR(particle.copy());
    }


    template<auto type, typename Particles, typename Electromag>
    void move_aos_pc_cpu(Particles& particles, Electromag const& em) const
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;

        auto per_particle = [&]<bool accelerate = false>() mutable {
            for (auto const& bix : particles.local_box())
                for (std::size_t idx = 0; idx < particles.size(bix); ++idx)
                {
                    auto& particle = particles(bix)[idx];
                    if constexpr (accelerate)
                        boris_accelerate(particle, interpolator_.m2p(particle, em, layout_), dto2m);
                    auto const& newCell = advancePosition_<alloc_mode>(particle, halfDtOverDl_);
                    if (newCell != particle.iCell())
                    {
                        particles.icell_changer(particle, bix, idx, newCell);
                        particle.iCell() = newCell;
                    }
                }
        };

        PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu");

        {
            PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu::advance_0");
            per_particle();
        }
        ParticleArrayService::sync<0, type>(particles);

        {
            PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu::advance_accelerate");
            per_particle.template operator()<true>();
        }
        ParticleArrayService::sync<1, type>(particles);
    }


    template<auto type, typename Particles, typename Electromag>
    void move_aos_pc_gpu(Particles& particles, Electromag const& em) const
    {
        PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu");

        if constexpr (any_in(Particles::impl_v, 0, 1))
            move_aos_pc_gpu_impl<type>(particles, em);

        else if constexpr (any_in(Particles::impl_v, 2))
            move_aos_pc_gpu_impl_2<type>(particles, em);

        else
            throw std::runtime_error("Unhandled AoS PC Boris impl");
    }

    template<auto type, typename Particles, typename Electromag>
    void move_aos_pc_gpu_impl(Particles& particles, Electromag const& em) const
    {
        PHARE_WITH_MKN_GPU({
            auto constexpr static alloc_mode = Particles::alloc_mode;
            auto const dto2m_                = dto2m;
            auto const layout                = layout_;
            auto const halfDtOverDl          = halfDtOverDl_;
            auto view                        = *particles;

            auto per_particle = [=] _PHARE_ALL_FN_<bool accelerate = false>() mutable
            {
                Interpolator interp;
                auto it        = view[mkn::gpu::idx()];
                auto& particle = *it;
                if constexpr (accelerate)
                    boris_accelerate(particle, interp.m2p(particle, em, layout), dto2m_);
                auto const& newCell = advancePosition_<alloc_mode>(particle, halfDtOverDl);
                if (!array_equals(newCell, particle.iCell()))
                {
                    it.icell_changer(newCell);
                    particle.iCell() = newCell;
                }
            };

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_0/1");
                mkn::gpu::GDLauncher{particles.size()}(per_particle);
            }

            if constexpr (any_in(Particles::impl_v, 1)) // on gpu particle swap/sync
                ParticleArrayService::sync<0, type>(view);
            ParticleArrayService::sync<0, type>(particles);

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_0/1");
                mkn::gpu::GDLauncher{particles.size()}(
                    [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });
            }

            if constexpr (any_in(Particles::impl_v, 1))
                ParticleArrayService::sync<1, type>(view);
            ParticleArrayService::sync<1, type>(particles);
        })
    }

    template<auto type, typename Particles, typename Electromag>
    void move_aos_pc_gpu_impl_2(Particles& particles, Electromag const& em) const
    {
        PHARE_WITH_MKN_GPU({
            auto constexpr static alloc_mode = Particles::alloc_mode;
            using lobox_t                    = Particles::lobox_t;
            using Launcher                   = gpu::BoxCellNLauncher<lobox_t>;

            auto const dto2m_       = dto2m;
            auto const layout       = layout_;
            auto const halfDtOverDl = halfDtOverDl_;
            auto view               = *particles;
            auto const& lobox       = view.local_box();

            auto per_particle = [=] _PHARE_ALL_FN_<bool accelerate = false>() mutable
            {
                auto const blockidx  = Launcher::block_idx();
                auto const threadIdx = Launcher::thread_idx();
                auto const& locell   = *(lobox.begin() + blockidx);
                auto& parts          = view(locell);
                if (threadIdx >= parts.size())
                    return;
                auto& particle = parts[threadIdx];
                Interpolator interp;
                if constexpr (accelerate)
                    boris_accelerate(particle, interp.m2p(particle, em, layout), dto2m_);
                auto const& newCell = advancePosition_<alloc_mode>(particle, halfDtOverDl);
                if (!array_equals(newCell, particle.iCell()))
                {
                    view.icell_changer(particle, locell, threadIdx, newCell);
                    particle.iCell() = newCell;
                }
            };

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_2");
                auto lobox = view.local_box();
                Launcher{lobox, particles.max_size()}(per_particle);
            }
            ParticleArrayService::sync<0, type>(view);      // no realloc
            ParticleArrayService::sync<0, type>(particles); // allows relloc

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_2");
                auto lobox = view.local_box();
                Launcher{lobox, particles.max_size()}(
                    [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });
            }
            ParticleArrayService::sync<1, type>(view);      // no realloc
            ParticleArrayService::sync<1, type>(particles); // allows relloc
        })
    }



    mutable Interpolator interpolator_;
    GridLayout const& layout_;
    std::array<double, dim> const halfDtOverDl_;
    double const mass_ = 0;
    double const dt_   = 0;
    double const dto2m = 0.5 * dt_ / mass_;
};

} // namespace PHARE::core


#endif /* PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP */
