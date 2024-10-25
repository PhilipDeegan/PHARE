#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP


#include "core/data/particles/particle_array_def.hpp"
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
    template<auto type, typename ModelViews>
    void particleToMesh(ModelViews& views, double coef = 1.) _PHARE_ALL_FN_
    {
        //
    }

    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    void particleToMesh(Particles const& particles, GridLayout const& layout, Field& density,
                        VecField& flux, double coef = 1.) _PHARE_ALL_FN_
    {
        static_assert(Particles::storage_mode == StorageMode::SPAN);
        auto constexpr static alloc_mode = ParticleArray_t::alloc_mode;
        auto constexpr static impl       = ParticleArray_t::impl;

        PHARE_LOG_SCOPE(1, "Interpolating::particleToMesh");

        if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        {
            if constexpr (alloc_mode == AllocatorMode::CPU)
            {
                for (auto const& bix : particles.local_box())
                    for (auto const& p : particles(bix))
                        interp_.particleToMesh(p, density, flux, layout, coef);
            }
            else if (alloc_mode == AllocatorMode::GPU_UNIFIED and impl < 2)
            {
                static_assert(atomic_ops, "GPU must be atomic");
                PHARE_WITH_MKN_GPU(
                    mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
                        Interpolator_t{}.particleToMesh(particles[mkn::gpu::idx()], density, flux,
                                                        layout, coef);
                    }); //
                )
            }
            else if (alloc_mode == AllocatorMode::GPU_UNIFIED and impl == 2)
            {
                PHARE_WITH_MKN_GPU({ //
                    using lobox_t  = Particles::lobox_t;
                    using Launcher = gpu::BoxCellNLauncher<lobox_t>;
                    auto lobox     = particles.local_box();
                    Launcher{lobox, particles.max_size()}([=] _PHARE_ALL_FN_() mutable {
                        box_kernel(particles, layout, flux, density);
                    });
                })
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
            static_assert(atomic_ops, "GPU must be atomic");
            PHARE_WITH_MKN_GPU( //
                mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
                    auto particle = particles.begin() + mkn::gpu::idx();
                    Interpolator_t{}.particleToMesh(deref(particle), density, flux, layout, coef);
                }); //
            )
        }
        else
            throw std::runtime_error("fail");
    }

    template<bool in_box = false, typename Particles, typename GridLayout, typename VecField,
             typename Field>
    static void box_kernel(Particles& particles, GridLayout& layout, VecField& flux, Field& density,
                           double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU
        using lobox_t         = Particles::lobox_t;
        using Launcher        = gpu::BoxCellNLauncher<lobox_t>;
        auto const& lobox     = particles.local_box();
        auto const& blockidx  = Launcher::block_idx();
        auto const& threadIdx = Launcher::thread_idx();
        auto& parts           = particles(*(lobox.begin() + blockidx));

        auto const interp = [&]() {
            Interpolator_t{}.particleToMesh(parts[threadIdx], density, flux, layout, coef);
        };

        if (threadIdx < parts.size())
        {
            if constexpr (!in_box)
            {
                interp();
            }
            else
            {
                if (isIn(parts[threadIdx], particles.box()))
                    interp();
            }
        }
#endif // PHARE_HAVE_MKN_GPU
    }


    Interpolator_t interp_;
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
