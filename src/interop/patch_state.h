
#ifndef PHARE_INTEROP_PATCH_STATE_H
#define PHARE_INTEROP_PATCH_STATE_H

#include <vector>
#include "phare_core.h"


namespace PHARE
{
template<typename GridLayout_>
struct PatchState
{
    using Float                     = double;
    using GridLayout                = GridLayout_;
    static constexpr auto dimension = GridLayout::dimension;
    using HybridQuantity            = core::HybridQuantity;
    using ParticleArray_t =
        typename core::PHARE_Types<dimension, GridLayout::interp_order>::ParticleArray_t;

    template<typename State>
    PatchState(GridLayout const& gridLayout, State& state)
        : layout{gridLayout}
    {
        auto& E = state.electromag.E;
        electromag.emplace_back(view(E[0].data(), HybridQuantity::Scalar::Ex));
        electromag.emplace_back(view(E[1].data(), HybridQuantity::Scalar::Ey));
        electromag.emplace_back(view(E[2].data(), HybridQuantity::Scalar::Ez));

        auto& B = state.electromag.B;
        electromag.emplace_back(view(B[0].data(), HybridQuantity::Scalar::Bx));
        electromag.emplace_back(view(B[1].data(), HybridQuantity::Scalar::By));
        electromag.emplace_back(view(B[2].data(), HybridQuantity::Scalar::Bz));

        for (auto& pop : state.ions)
        {
            ions.emplace_back(&pop.domainParticles());
            assert(ions.back()->size()); // ?
            masses.emplace_back(pop.mass());

            density.emplace_back(view(pop.density().data(), HybridQuantity::Scalar::rho));
            assert(density.back().data() == pop.density().data());

            auto& F = pop.flux();
            flux.emplace_back(view(F[0].data(), HybridQuantity::Scalar::rho));
            flux.emplace_back(view(F[1].data(), HybridQuantity::Scalar::rho));
            flux.emplace_back(view(F[2].data(), HybridQuantity::Scalar::rho));
        }
    }

    template<typename T, typename Qty>
    static auto view(GridLayout_ const& layout, T* ptr, Qty qty)
    {
        return core::NdArrayView<dimension, T>{ptr, layout.allocSize(qty)};
    }
    template<typename T, typename Qty>
    auto view(T* ptr, Qty qty)
    {
        return view(layout, ptr, qty);
    }


    GridLayout const layout;
    std::vector<core::NdArrayView<dimension, Float>> electromag;
    std::vector<core::NdArrayView<dimension, Float>> density, flux;
    std::vector<ParticleArray_t*> ions;
    std::vector<double> masses;
};

} // namespace PHARE


#if defined(HAVE_UMPIRE)
#include "kul/gpu.hpp"
namespace PHARE
{
template<typename PatchState0, typename PatchState1>
void load_to_umpire_device_mem(PatchState0 const&& src, PatchState1&& dst)
{
    auto send = [](auto& src_vec, auto& dst_vec) {
        for (std::size_t i = 0; i < src_vec.size(); ++i)
            KUL_GPU_NS::send(dst_vec[i].data(), src_vec[i].data(), src_vec[i].size());
    };

    send(src.electromag, dst.electromag);
    send(src.flux, dst.flux);
    send(src.ions, dst.ions);
}

template<typename PatchState0, typename PatchState1>
void load_from_umpire_device_mem(PatchState0 const&& src, PatchState1&& dst)
{
    auto take = [](auto& src_vec, auto& dst_vec) {
        for (std::size_t i = 0; i < src_vec.size(); ++i)
            KUL_GPU_NS::take(src_vec[i].data(), dst_vec[i].data(), src_vec[i].size());
    };

    take(src.electromag, dst.electromag);
    take(src.flux, dst.flux);
    take(src.ions, dst.ions);
}

} // namespace PHARE

#endif /* HAVE_UMPIRE */
#endif /* PHARE_INTEROP_PATCH_STATE_H */
