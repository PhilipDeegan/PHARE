#ifndef PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP
#define PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP

// #include "core/utilities/kernels.hpp"
#include "mkn_boris.hpp"

namespace PHARE::core::boris
{




template<typename... Args>
auto mover(Args&&... args)
{
    if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(return MultiBoris{args...});
    }
    else
        throw std::runtime_error("Vector::copy NO ALTERNATIVE");
}


#if 0
auto per_particle = [=] _PHARE_DEV_FN_(auto const& i) mutable {
    using ParticleArray_v            = typename MultiBoris<ModelViews>::ParticleArray_v;
    using BoxLauncher    = gpu::BoxCellNLauncher<Box<std::uint32_t, dim>, false>;
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

    if constexpr (any_in(ParticleArray_v::layout_mode, LayoutMode::SoA,
    LayoutMode::SoAPC))
    {
        detail::SoAZipParticle particle{parts, threadIdx};
        auto const og_iCell = particle.iCell();

        particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
        {
            Interpolator interp;
            boris::accelerate(particle, interp.m2p(particle, emps[i], layout), dto2m);
        }

        particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

        if (!array_equals(particle.iCell(), og_iCell))
            view.icell_changer(particle, locell, threadIdx, particle.iCell());
    }
    else
    {
        auto& particle      = parts[threadIdx];
        auto const og_iCell = particle.iCell();

        particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
        {
            Interpolator interp;
            boris::accelerate(particle, interp.m2p(particle, emps[i], layout), dto2m);
        }

        particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

        if (!array_equals(particle.iCell(), og_iCell))
            view.icell_changer(particle, locell, threadIdx, particle.iCell());
    }
};


// OR

// std::size_t const ppc = pps[0]({2, 2, 2}).size();
auto per_chunk = [=] _PHARE_DEV_FN_(auto const& i) mutable {
    using ParticleArray_v            = typename MultiBoris<ModelViews>::ParticleArray_v;
    // auto const tidx = mkn::gpu::idx();
    // assert(tidx < 65);

    // auto const blockidx
    //     = (8 * 8 * ((tidx / 16) + 2)) + 16 + (((tidx % 16) / 4) * 8) + 2 + (tidx % 4);
    auto const t_x = threadIdx.x + 2;
    auto const t_y = threadIdx.y + 2;
    auto const t_z = threadIdx.z + 2;

    auto const& dto2m  = dto2mspp[i];
    auto const& layout = layoutps[i];
    auto& view         = pps[i];
    // auto const& lobox  = view.local_box();
    // auto const& locell = *(lobox.begin() + blockidx);
    std::array<std::uint32_t, 3> locell{t_x, t_y, t_z};
    auto& parts = *view().data();

    for (std::size_t pidx = 0; pidx < ppc; ++pidx)
    {
        if constexpr (any_in(ParticleArray_v::layout_mode, LayoutMode::SoA,
                             LayoutMode::SoAPC))
        {
            detail::SoAZipParticle particle{parts, pidx};
            auto const og_iCell = particle.iCell();

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
            {
                Interpolator interp;
                boris::accelerate(particle, interp.m2p(particle, emps[i], layout),
                dto2m);
            }

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

            if (!array_equals(particle.iCell(), og_iCell))
                view.icell_changer(particle, locell, pidx, particle.iCell());
        }
        else
        {
            auto& particle      = parts[pidx];
            auto const og_iCell = particle.iCell();

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);
            {
                Interpolator interp;
                boris::accelerate(particle, interp.m2p(particle, emps[i], layout),
                dto2m);
            }

            particle.iCell() = boris::advance<alloc_mode>(particle, halfDtOverDl[i]);

            if (!array_equals(particle.iCell(), og_iCell))
                view.icell_changer(particle, locell, pidx, particle.iCell());
        }
    }
};
#endif // 0


} // namespace PHARE::core::boris


#endif /*PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP*/
