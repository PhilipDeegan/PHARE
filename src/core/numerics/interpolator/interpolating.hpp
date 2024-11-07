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
                        auto it = particles[mkn::gpu::idx()];
                        Interpolator_t{}.particleToMesh(*it, density, flux, layout, coef);
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

    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel(Particles& particles, GridLayout& layout, VecField& flux,
                             Field& density, double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU
        auto constexpr static nghosts = GridLayout::nbrGhosts();
        extern __shared__ double data[];
        auto const& lobox = particles.local_box();
        auto const ziz    = 9 * 9 * 9;

        auto const t_x = threadIdx.x;
        auto const t_y = threadIdx.y;
        auto const t_z = threadIdx.z;

        core::Point<std::uint32_t, 3> const tcell{t_x, t_y, t_z};
        core::Point<std::uint32_t, 3> locell = tcell + nghosts;

        auto rho = make_array_view(&data[ziz * 0], density.shape());
        auto fx  = make_array_view(&data[ziz * 1], flux[0].shape());
        auto fy  = make_array_view(&data[ziz * 2], flux[1].shape());
        auto fz  = make_array_view(&data[ziz * 3], flux[2].shape());

        // global mem pointer backps
        auto const r0 = density.data();
        auto const f0 = flux[0].data();
        auto const f1 = flux[1].data();
        auto const f2 = flux[2].data();

        density.setBuffer(rho.data());
        flux[0].setBuffer(fx.data());
        flux[1].setBuffer(fy.data());
        flux[2].setBuffer(fz.data());

        for (std::size_t i = 0; i < 9 * 9 * 9; ++i)
        {
            auto const cell = *(lobox.begin() + i);
            density(cell)   = 0;
            flux[0](cell)   = 0;
            flux[1](cell)   = 0;
            flux[2](cell)   = 0;
        }

        __syncthreads();

        Interpolator_t interp;
        auto& parts = particles(locell);
        for (std::size_t i = 0; i < parts.size(); ++i)
            interp.particleToMesh(parts[i], density, flux, layout, coef);

        __syncthreads();

        density.setBuffer(r0);
        flux[0].setBuffer(f0);
        flux[1].setBuffer(f1);
        flux[2].setBuffer(f2);

        for (std::size_t i = 0; i < 9 * 9 * 9; ++i)
        {
            auto const cell = *(lobox.begin() + i);
            density(cell)   = rho(cell);
            flux[0](cell)   = fx(cell);
            flux[1](cell)   = fy(cell);
            flux[2](cell)   = fz(cell);
        }

#endif // PHARE_HAVE_MKN_GPU
    }


    Interpolator_t interp_;
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
