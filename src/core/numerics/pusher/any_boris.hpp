#ifndef PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP
#define PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP

#include "core/utilities/kernels.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"


#include <cmath>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>


namespace PHARE::core
{



template<std::size_t dim, typename Interpolator, typename GridLayout>
class AnyBorisPusher
{
    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }

public:
    AnyBorisPusher(                        //
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
            auto view = *particles;
            PHARE_ASSERT(view.size());
            kernel::launch(particles.size(), [=, layout = layout_, halfDtOverDl = halfDtOverDl_,
                                              dto2m_ = dto2m] _PHARE_ALL_FN_() mutable {
                Interpolator interp;
                auto particle = view.begin() + kernel::idx();
                pp<alloc_mode>(particle, em, interp, layout, halfDtOverDl, dto2m_);
            });
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
                auto it        = view[kernel::idx()];
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
                kernel::launch(particles.size(), per_particle);
                // mkn::gpu::GDLauncher{particles.size()}(per_particle);
            }

            if constexpr (any_in(Particles::impl_v, 1)) // on gpu particle swap/sync
                ParticleArrayService::sync<0, type>(view);
            ParticleArrayService::sync<0, type>(particles);

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_0/1");
                kernel::launch(particles.size(), [=] _PHARE_ALL_FN_() mutable {
                    per_particle.template operator()<true>();
                });
                // mkn::gpu::GDLauncher{particles.size()}(
                //     [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });
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
            // using Launcher                   = gpu::BoxCellNLauncher<lobox_t>;

            auto const dto2m_       = dto2m;
            auto const layout       = layout_;
            auto const halfDtOverDl = halfDtOverDl_;
            auto view               = *particles;
            auto const& lobox       = view.local_box();

            auto per_particle = [=] _PHARE_ALL_FN_<bool accelerate = false>() mutable
            {
                auto const blockidx  = kernel::block_idx();
                auto const threadIdx = kernel::thread_idx();
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
                // auto lobox = view.local_box();

                kernel::launch(view.local_box(), particles.max_size(), per_particle);

                // Launcher{lobox, particles.max_size()}(per_particle);
            }
            ParticleArrayService::sync<0, type>(view);      // no realloc
            ParticleArrayService::sync<0, type>(particles); // allows relloc

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_2");
                // auto lobox = view.local_box();
                kernel::launch(
                    view.local_box(), particles.max_size(),
                    [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });

                // Launcher{lobox, particles.max_size()}(
                //     [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });
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
