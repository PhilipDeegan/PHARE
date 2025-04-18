#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/utilities/kernels/launchers.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "interpolator.hpp"


namespace PHARE::core
{

// simple facade to launch e.g. GPU kernels if needed
//  and to not need to modify the interpolator much for gpu specifically

template<typename ParticleArray_t, std::size_t interpOrder, bool atomic_ops,
         typename Interpolator_t
         = Interpolator<ParticleArray_t::dimension, interpOrder, atomic_ops>>
class Interpolating
{
    auto constexpr static dim = ParticleArray_t::dimension;

public:
    template<typename Particles, typename VecField, typename GridLayout, typename Field>
    inline void operator()(Particles const& particles, Field& density, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
    {
        using enum LayoutMode;
        if constexpr (any_in(ParticleArray_t::layout_mode, AoSMapped))
            Interpolator_t{}(particles, density, flux, layout, coef);
        else
            particleToMesh(particles, layout, density, flux, coef);
    }


    template<typename GridLayout, typename VecField, typename Field>
    void particleToMesh(ParticleArray_t const& particles, GridLayout const& layout, Field& density,
                        VecField& flux, double coef = 1.) _PHARE_ALL_FN_
    {
        using Particles = ParticleArray_t;
        // static_assert(Particles::storage_mode == StorageMode::SPAN);
        auto constexpr static alloc_mode = ParticleArray_t::alloc_mode;
        auto constexpr static impl       = ParticleArray_t::impl;

        PHARE_LOG_SCOPE(1, "Interpolating::particleToMesh");

        using enum LayoutMode;
        if constexpr (any_in(Particles::layout_mode, AoSTS))
        {
            if constexpr (alloc_mode == AllocatorMode::CPU)
            {
                for (auto const& tile : particles())
                    for (auto const& p : tile())
                        interp_.particleToMesh(p, density, flux, layout, coef);
            }
            else
            {
                throw std::runtime_error("fail");
            }
        }
        else if constexpr (Particles::layout_mode == AoSPC)
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
            if constexpr (any_in(Particles::layout_mode, SoATS))
                for (auto const& tile : particles())
                    for (std::size_t i = 0; i < tile().size(); ++i)
                        interp_.particleToMesh(tile()[i], density, flux, layout, coef);
            else if constexpr (any_in(Particles::layout_mode, SoAVX))
            {
                for (auto const& it : particles)
                    interp_.particleToMesh(it, density, flux, layout, coef);
            }
            else
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
#if PHARE_HAVE_MKN_GPU_HW
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
#endif // PHARE_HAVE_MKN_GPU_HW
    }

    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel_ts(Particles& particles, GridLayout& layout, VecField& flux,
                                Field& density, double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW

        Point<std::uint32_t, 3> const rcell{threadIdx.z, threadIdx.y, threadIdx.x};
        extern __shared__ double data[];

        auto& tile        = *particles().at(rcell + 2);
        auto& parts       = tile();
        auto const each   = parts.size() / 8;
        auto const tile_p = rcell % 2;
        auto const pidx   = tile_p[2] + tile_p[1] * 2 + tile_p[0] * 2 * 2;
        auto const ji     = [&](auto const i) { return i * 8 + pidx; };
        auto doX          = [&]<std::uint8_t IDX = 0>(auto& feeld) {
            auto v = make_array_view(&data[0], feeld.shape());
            if (mkn::gpu::idx() == 0)
                for (auto& e : v)
                    e = 0;
            Interpolator_t interp;
            std::size_t i = 0, j = ji(i);
            __syncthreads();
            for (; i < each; ++i, j = ji(i))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            if (pidx < parts.size() - (8 * each))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            __syncthreads();
            if (mkn::gpu::idx() == 0)
                for (i = 0; i < v.size(); ++i)
                    feeld.data()[i] = v.data()[i];
        };

        doX(density);
        doX.template operator()<1>(flux[0]);
        doX.template operator()<2>(flux[1]);
        doX.template operator()<3>(flux[2]);


#endif // PHARE_HAVE_MKN_GPU_HW
    }

    template<typename Field, typename Particles, typename GridLayout>
    static void on_tiles(Particles& particles, GridLayout& layout, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW

        static_assert(atomic_ops, "GPU must be atomic");
        using TensorField_t = basic::TensorField<typename Field::Super::Super, 1>;
        extern __shared__ double data[];
        // printf("L:%d i %llu \n", __LINE__, blockIdx.x);
        auto const ziz  = 9 * 9 * 9;
        auto& tile      = particles()[blockIdx.x];
        auto& parts     = tile();
        auto const tbox = tile.field_box(); //(*tile).unsafe_intersection(particles.box());
        auto const tidx = threadIdx.x;
        auto const ws   = kernel::warp_size();
        auto const each = parts.size() / ws;
        auto const& [density, flux0, flux1, flux2] = tile.fields();
        auto const r0                              = density;

        TensorField_t flux{flux0, flux1, flux2};

        auto rho = make_array_view(&data[ziz * 0], density.shape());
        auto fx  = make_array_view(&data[ziz * 1], flux[0].shape());
        auto fy  = make_array_view(&data[ziz * 2], flux[1].shape());
        auto fz  = make_array_view(&data[ziz * 3], flux[2].shape());

        density.reset(rho);
        flux[0].reset(fx);
        flux[1].reset(fy);
        flux[2].reset(fz);

        if (tidx == 0)
            for (std::size_t i = 0; i < ziz * 4; ++i)
                data[i] = 0;

        Interpolator_t interp;
        std::size_t pid = 0;

        __syncthreads();
        for (; pid < each; ++pid)
            interp.particleToMesh(parts[(pid * ws) + tidx], density, flux, layout, tbox, coef);
        if (tidx < parts.size() - (ws * each))
            interp.particleToMesh(parts[(pid * ws) + tidx], density, flux, layout, tbox, coef);
        __syncthreads();

        density.reset(r0);

        if (tidx == 0)
        {
            for (std::size_t i = 0; i < rho.size(); ++i)
            {
                density.data()[i] = rho.data()[i];
                flux0.data()[i]   = fx.data()[i];
                flux1.data()[i]   = fy.data()[i];
                flux2.data()[i]   = fz.data()[i];
            }
        };

#endif
    }


    template<typename Field, typename Particles, typename GridLayout>
    static void on_cpu_tiles(Particles& particles, GridLayout& layout,
                             double coef = 1.) _PHARE_ALL_FN_
    {
        using TensorField_t = basic::TensorField<typename Field::Super, 1>;

        for (auto& tile : particles())
        {
            auto& parts                                = tile();
            auto const tbox                            = tile.field_box();
            auto const& [density, flux0, flux1, flux2] = tile.fields();
            auto const r0                              = density;

            TensorField_t flux{flux0, flux1, flux2};
            Interpolator_t interp;
            std::size_t pid = 0;

            interp.particleToMesh(parts[pid], density, flux, layout, tbox, coef);
        }
    }



    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel_ts_all(Particles& particles, GridLayout& layout, VecField& flux,
                                    Field& density, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW

        Point<std::uint32_t, 3> const rcell{threadIdx.z, threadIdx.y, threadIdx.x};
        extern __shared__ double data[];
        auto const ziz    = 9 * 9 * 9;
        auto& tile        = *particles().at(rcell + 2);
        auto& parts       = tile();
        auto const each   = parts.size() / 8;
        auto const tile_p = rcell % 2;
        auto const pidx   = tile_p[2] + tile_p[1] * 2 + tile_p[0] * 2 * 2;
        auto const ji     = [&](auto const i) { return i * 8 + pidx; };
        auto doX          = [](auto const& v, auto& feeld) {
            if (mkn::gpu::idx() == 0)
                for (std::size_t i = 0; i < v.size(); ++i)
                    feeld.data()[i] = v.data()[i];
        };

        auto const r0 = density.data();
        auto const f0 = flux[0].data();
        auto const f1 = flux[1].data();
        auto const f2 = flux[2].data();

        auto rho = make_array_view(&data[ziz * 0], density.shape());
        auto fx  = make_array_view(&data[ziz * 1], flux[0].shape());
        auto fy  = make_array_view(&data[ziz * 2], flux[1].shape());
        auto fz  = make_array_view(&data[ziz * 3], flux[2].shape());

        density.setData(rho.data());
        flux[0].setData(fx.data());
        flux[1].setData(fy.data());
        flux[2].setData(fz.data());

        auto const z0 = [](auto& arr) {
            for (auto& e : arr)
                e = 0;
        };
        if (mkn::gpu::idx() == 0)
            [&](auto&... arr) { (z0(arr), ...); }(rho, fx, fy, fz);

        Interpolator_t interp;
        std::size_t i = 0, j = ji(i);

        __syncthreads();
        for (; i < each; ++i, j = ji(i))
            interp.particleToMesh(parts[j], density, flux, layout, coef);
        if (pidx < parts.size() - (8 * each))
            interp.particleToMesh(parts[j], density, flux, layout, coef);
        __syncthreads();

        density.setData(r0);
        flux[0].setData(f0);
        flux[1].setData(f1);
        flux[2].setData(f2);


        doX(rho, density);
        doX(fx, flux[0]);
        doX(fy, flux[1]);
        doX(fz, flux[2]);


#endif // PHARE_HAVE_MKN_GPU_HW
    }


    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel(Particles& particles, GridLayout& layout, VecField& flux,
                             Field& density, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        auto constexpr static nghosts = GridLayout::nbrGhosts();
        extern __shared__ double data[];
        auto const& lobox = particles.local_box();
        auto const ziz    = 9 * 9 * 9;

        auto const t_x = threadIdx.x;
        auto const t_y = threadIdx.y;
        auto const t_z = threadIdx.z;

        core::Point<std::uint32_t, 3> const tcell{t_x, t_y, t_z};
        core::Point<std::uint32_t, 3> const locell = tcell + nghosts;

        {
            Interpolator_t interp;
            auto& parts = particles(locell);

            { // rho
                auto const r0 = density.data();
                auto v        = make_array_view(&data[0], density.shape());
                density.setData(v.data());
                if (mkn::gpu::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], v);
                }
                __syncthreads();
                density.setData(r0);
                if (mkn::gpu::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        density.data()[i] = v.data()[i];
            }

            {
                auto const f0 = flux[0].data();
                auto v        = make_array_view(&data[0], flux[0].shape());
                flux[0].setData(v.data());
                if (mkn::gpu::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[0]);
                }
                __syncthreads();
                flux[0].setData(f0);
                if (mkn::gpu::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[0].data()[i] = v.data()[i];
            }

            {
                auto const f1 = flux[1].data();
                auto v        = make_array_view(&data[0], flux[1].shape());
                flux[1].setData(v.data());
                if (mkn::gpu::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[1]);
                }
                __syncthreads();
                flux[1].setData(f1);
                if (mkn::gpu::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[1].data()[i] = v.data()[i];
            }

            {
                auto const f2 = flux[2].data();
                auto v        = make_array_view(&data[0], flux[2].shape());
                flux[2].setData(v.data());
                if (mkn::gpu::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[2]);
                }
                __syncthreads();
                flux[2].setData(f2);
                if (mkn::gpu::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[2].data()[i] = v.data()[i];
            }
        }


#endif // PHARE_HAVE_MKN_GPU_HW
    }

    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void ts_reducer(Particles& particles, GridLayout const& layout, VecField& flux,
                           Field& density) _PHARE_ALL_FN_
    {
        PHARE_LOG_LINE_SS("");
        auto const safe_box = particles.safe_box();

        for (auto& tile : particles())
        {
            auto const& tile_box         = tile.field_box();
            auto const& [r0, f0, f1, f2] = tile.fields();
            for (auto const& bix : tile_box)
            {
                auto const plix = (bix - safe_box.lower()).as_unsigned();
                auto const tlix = (bix - tile_box.lower()).as_unsigned();

                density(plix) += r0(tlix);
                flux[0](plix) += f0(tlix);
                flux[1](plix) += f1(tlix);
                flux[2](plix) += f2(tlix);
            }
        }
    }


    Interpolator_t interp_;
};

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
struct MomentumTensorInterpolating
{
    using Interpolator_t = MomentumTensorInterpolator<dim, interpOrder, atomic_ops>;

    Interpolator_t interp_;

public:
    template<typename Particles_t, typename TensorField, typename GridLayout>
    inline void operator()(Particles_t& particles, TensorField& momentumTensor,
                           GridLayout const& layout, double mass = 1.)
    {
        using enum LayoutMode;

        if constexpr (any_in(Particles_t::layout_mode, AoSMapped))
            interp_(particles, momentumTensor, layout, mass);
        else if constexpr (any_in(Particles_t::layout_mode, AoSTS))
        {
            for (auto const& tile : particles())
                interp_(tile(), momentumTensor, layout, mass);
        }
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
