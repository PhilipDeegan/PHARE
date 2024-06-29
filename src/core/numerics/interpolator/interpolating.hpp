#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "interpolator.hpp"



namespace PHARE::core
{

// simple facade to launch e.g. GPU kernels if needed
//  and to not need to modify the interpolator much for gpu specifically

template<typename ParticleArray_t, std::size_t interpOrder, bool atomic_ops>
class Interpolating
{
    auto constexpr static dim = ParticleArray_t::dimension;
    using Interpolator_t      = Interpolator<dim, interpOrder, atomic_ops>;

public:
    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    void particleToMesh(Particles const& particles, GridLayout const& layout, Field& density,
                        VecField& flux, double coef = 1.) _PHARE_ALL_FN_
    {
        static_assert(Particles::storage_mode == StorageMode::SPAN);
        auto constexpr static alloc_mode = ParticleArray_t::alloc_mode;

        if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        {
            if constexpr (alloc_mode == AllocatorMode::CPU)
            {
                // interp_(particles, density, flux, layout, coef);
                // Interpolator<dim, interpOrder, false> _interp; // force CPU
                for (auto const& bix : particles.local_box())
                    for (auto const& p : particles(bix))
                        interp_.particleToMesh(p, density, flux, layout, coef);
            }
            else if (alloc_mode == AllocatorMode::GPU_UNIFIED)
            {
                particles.print();
                particles.check();
                assert(mkn::gpu::Pointer{density.data()}.is_device_ptr());
                assert(mkn::gpu::Pointer{flux[0].data()}.is_device_ptr());
                assert(mkn::gpu::Pointer{&particles[0].delta()[0]}.is_device_ptr());
                static_assert(atomic_ops, "GPU must be atomic");
                PHARE_WITH_MKN_GPU(
                    mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
                        Interpolator_t{}.particleToMesh(particles[mkn::gpu::idx()], density, flux,
                                                        layout, coef);
                    }); //
                )
            }
            else
                throw std::runtime_error("fail");
        }
        else if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            interp_(particles, density, flux, layout, coef);
        }
        else if (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            static_assert(atomic_ops, "GPU must be atomic"); //
            // auto& rho_view  = *density;
            // auto& flux_view = *flux;
            // auto ps = particles.view();
            PHARE_WITH_MKN_GPU( //
                mkn::gpu::GDLauncher{particles.size()}(
                    [=, interp = interp_] _PHARE_ALL_FN_() mutable {
                        auto i        = interp; // needs own state
                        auto particle = particles.begin() + mkn::gpu::idx();
                        interp.particleToMesh(deref(particle), density, flux, layout, coef);
                    }); //
            )
        }
        else
            throw std::runtime_error("fail");
    }


    // void m2p(){};

    Interpolator_t interp_;
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
