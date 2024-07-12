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


#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"

#include "core/numerics/pusher/pusher.hpp"

#include "core/numerics/pusher/boris_simpler.hpp"


namespace PHARE::core
{
template<std::size_t dim, typename Interpolator, typename GridLayout>
class MultiBorisPusher
{
    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }

public:
    MultiBorisPusher(double dt) : dt_{dt}
    {
    }

    template<auto alloc_mode, auto layout_mode, typename ModelViews>
    void operator()(ModelViews & views) const
    {
        move<alloc_mode, layout_mode, ParticleType::Domain>(views);
    }

    template<auto alloc_mode, auto layout_mode, auto type, typename ModelViews>
    void push(ModelViews & views) const
    {
        move<alloc_mode, layout_mode, type>(views);
    }


    template<auto alloc_mode, auto layout_mode, auto type, typename ModelViews>
    void move(ModelViews & views) const _PHARE_ALL_FN_
    {
        static_assert(Particles::layout_mode == LayoutMode::AoSPC && alloc_mode == AllocatorMode::GPU_UNIFIED);

        // PHARE_LOG_LINE_STR(particles.size());
        if (particles.size() == 0)
            return; // noop

        PHARE_LOG_SCOPE(1, "multiboris::move");

        if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        {
            if constexpr (alloc_mode == AllocatorMode::CPU)
                move_aos_pc_cpu<type>(particles, em);
            else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
                move_aos_pc_gpu<type>(particles, em);
            else
                throw std::runtime_error("no");
        }
        // else if constexpr (alloc_mode == AllocatorMode::CPU)
        // {
        //     for (auto& p : particles)
        //         pp<alloc_mode>(p, em, interpolator_, layout_, halfDtOverDl_, dto2m);
        // }
        // else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        // {
        //     PHARE_WITH_MKN_GPU(            //
        //         auto view = *particles;    //
        //         PHARE_ASSERT(view.size()); //
        //         mkn::gpu::GDLauncher{particles.size()}([=, layout = layout_,
        //                                                 halfDtOverDl = halfDtOverDl_,
        //                                                 dto2m_ = dto2m] _PHARE_ALL_FN_() mutable {
        //             Interpolator interp;
        //             auto particle = view.begin() + mkn::gpu::idx();
        //             pp<alloc_mode>(particle, em, interp, layout, halfDtOverDl, dto2m_);
        //         }); //
        //     )
        // }
        else
            throw std::runtime_error("no");
    }




    template<auto alloc_mode, auto layout_mode, auto type, typename ModelViews>
    void move_aos_pc_gpu(ModelViews & views) const
    {
        bool constexpr static sync_on_kernel_finish = false;
        using ModelView = typename ModelViews::value_type;
        using ParticleArray_t = typenmae ModelView::ParticleArray_t;

        using Tup_t = std::tuple<ParticleArray_t*, mkn::gpu::GDLauncher<sync_on_kernel_finish>, mkn::gpu::Stream>;

        auto constexpr static alloc_mode = Particles::alloc_mode;
        PHARE_LOG_SCOPE(1, "multiboris::move_aos_pc_gpu");

        PHARE_WITH_MKN_GPU({

          std::vector<Tup_t> locals;

          for(auto & patch : views){
            for(auto& pop : *patch.ions){
              if constexpr(type == ParticleType::Domain){
                auto& tup = locals.emplace_back(&pop.domainParticles(), pop.domainParticles().size());
                std::get<1>(tup).s  = std::get<2>(tup)();
              }
            }
          }
        })

        // PHARE_WITH_MKN_GPU({
        //     PHARE_LOG_SCOPE(1, "multiboris::advance");
        //     // auto view = *particles;
        //     // view.cap(particles);

        //     mkn::gpu::GDLauncher{particles.size()}(
        //         [=, halfDtOverDl = halfDtOverDl_] _PHARE_ALL_FN_() mutable {
        //             auto particle       = view[mkn::gpu::idx()];
        //             auto const& newCell = advancePosition_<alloc_mode>(particle, halfDtOverDl);
        //             if (!array_equals(newCell, particle.iCell()))
        //             {
        //                 particle.icell_changer(newCell);
        //                 particle.iCell() = newCell;
        //             }
        //         });
        //     if constexpr (Particles::impl_v == 1)
        //         ParticleArrayService::sync<0, type>(view);
        // })

        // // ParticleArrayService::sync<0, type>(particles);

        // PHARE_WITH_MKN_GPU({
        //     PHARE_LOG_SCOPE(1, "multiboris::advance_accelerate");

        //     auto view = *particles;
        //     view.cap(particles);
        //     mkn::gpu::GDLauncher{particles.size()}([=, layout = layout_,
        //                                             halfDtOverDl = halfDtOverDl_,
        //                                             dto2m_       = dto2m] _PHARE_ALL_FN_() mutable {
        //         Interpolator interp;
        //         auto particle = view[mkn::gpu::idx()];
        //         boris_accelerate(particle, interp.m2p(particle, em, layout), dto2m_);
        //         auto const& newCell = advancePosition_<alloc_mode>(particle, halfDtOverDl);
        //         if (!array_equals(newCell, particle.iCell()))
        //         {
        //             particle.icell_changer(newCell);
        //             particle.iCell() = newCell;
        //         }
        //     });
        //     if constexpr (Particles::impl_v == 1)
        //         ParticleArrayService::sync<1, type>(view);
        // })

        // ParticleArrayService::sync<1, type>(particles);
    }


    mutable Interpolator interpolator_;
    GridLayout const& layout_;
    std::array<double, dim> const halfDtOverDl_;
    double const mass_ = 0;
    double const dt_   = 0;
    // double const dto2m = 0.5 * dt_ / mass_;
};

} // namespace PHARE::core


#endif /* PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP */
