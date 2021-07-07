#ifndef PHARE_GPU_HPP
#define PHARE_GPU_HPP

#include "kul/span.hpp"
#include "kul/gpu.hpp"
#include "kul/gpu/tuple.hpp"

#include "tests/simulator/per_test.h"

namespace PHARE::gpu
{
template<typename SimOpts>
struct PatchState
{
    using Float      = typename SimOpts::Float;
    using GridLayout = typename SimOpts::GridLayout_t;

    template<typename State>
    PatchState(GridLayout const& gridLayout, State const& state)
        : gridLayoutDAO{gridLayout}
    {
        for (auto const& pop : state.ions)
            ions.emplace_back(&pop.domainParticles());

        auto vecF = [&](auto const& EBxyz) {
            for (std::uint8_t i = 0; i < 3; i++)
                electromag.emplace_back(EBxyz[i].data(), EBxyz[i].size());
        };

        vecF(state.electromag.E);
        vecF(state.electromag.B);
    }

    typename GridLayout::Super gridLayoutDAO;
    std::vector<kul::Span<Float const, std::uint32_t>> electromag;
    std::vector<core::ParticleArray<Float, SimOpts::dimension>*> ions;
};


template<typename Float, std::uint8_t dim, bool GPU>
struct EMContiguousParticles : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = EMContiguousParticles<Float, dim, true>;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    EMContiguousParticles(std::uint32_t nbr) // bytes per particle
        : leaving{nbr}                       // 1
        , iCell{nbr * dim}                   // 4 * dim
        , delta{nbr * dim}                   // 4 * dim
        , weight{nbr}                        // 8
        , charge{nbr}                        // 8
        , v{nbr * 3}                         // 24
        , E{nbr * 3}                         // 24
        , B{nbr * 3}                         // 24
    {                                        // 1d  97 - 1e6 ~100MB / 1e7 ~1GB
    }                                        // 2d 105
                                             // 3d 113

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto add(core::ParticleArray<Float, dim> const& array, std::uint32_t i)
    {
        PHARE::core::EMContiguousParticles<Float, dim> particles{array};
        auto size      = array.size();
        auto dim_start = i * dim, _3_start = i * 3;
        auto dim_size = size * dim, _3_size = size * 3;

        iCell.send(&particles.iCell[0], dim_size, dim_start);
        delta.send(&particles.delta[0], dim_size, dim_start);
        weight.send(&particles.weight[0], size, i);
        charge.send(&particles.charge[0], size, i);
        v.send(&particles.v[0], _3_size, _3_start);
        E.send(&particles.E[0], _3_size, _3_start);
        B.send(&particles.B[0], _3_size, _3_start);

        return i + size;
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(leaving, iCell, delta, weight, charge, v, E, B);
    }

    container_t<bool> leaving;
    container_t<std::int32_t> iCell;
    container_t<Float> delta;
    container_t<Float> weight, charge, v, E, B;
};


template<typename Float, typename Span_ = kul::gpu::Span<Float const, std::uint32_t>>
struct FieldInterop
{
    using Span = Span_;

    FieldInterop() __device__ = default;

    auto& __device__ operator()(std::uint32_t i) { return ptrs[i]; }
    // auto operator()(std::uint32_t, std::uint32_t) { return 2; }
    // auto operator()(std::uint32_t, std::uint32_t, std::uint32_t) { return 3; }

    Span ptrs;
};

template<typename Float>
struct VecFieldInterop
{
    using Span = typename FieldInterop<Float>::Span;

    VecFieldInterop() __device__ = default;

    template<typename Electromags>
    static VecFieldInterop<Float> __device__ E(Electromags const& em)
    {
        return {{Span{em.Ex, em.info[0]}}, {Span{em.Ey, em.info[1]}}, {Span{em.Ez, em.info[2]}}};
    }
    template<typename Electromags>
    static VecFieldInterop<Float> __device__ B(Electromags const& em)
    {
        return {{Span{em.Bx, em.info[3]}}, {Span{em.By, em.info[4]}}, {Span{em.Bz, em.info[5]}}};
    }

    auto& __device__ getComponent(PHARE::core::Component XYZ)
    {
        if (XYZ == PHARE::core::Component::X)
            return x;
        else if (XYZ == PHARE::core::Component::Y)
            return y;
        return z;
    }

    FieldInterop<Float> x, y, z;
};


template<typename Float>
struct EMInterop
{
    template<typename Electromags>
    EMInterop(Electromags&& _em) __device__ : E{VecFieldInterop<Float>::E(_em)},
                                              B{VecFieldInterop<Float>::B(_em)}
    {
    }

    VecFieldInterop<Float> E, B;
};


template<typename SimOpts, bool GPU>
struct Electromags : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = Electromags<SimOpts, true>;
    using Float = typename SimOpts::Float;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static auto make_shared(std::vector<PatchState<SimOpts>> const& states)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());

        std::vector<std::uint32_t> emXYZ(6, 0), emInfo(n_states * 6);
        for (std::size_t j = 0; j < n_states; j++)
            for (std::size_t i = 0; i < states[j].electromag.size(); i++)
            {
                auto pos    = (j * 6) + i;
                emInfo[pos] = emXYZ[i];
                emXYZ[i] += states[j].electromag[i].size();
            }

        return std::make_shared<Electromags<SimOpts, GPU>>(states, emXYZ, emInfo);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    Electromags(std::vector<PatchState<SimOpts>> const& states, std::vector<std::uint32_t> const& v,
                std::vector<std::uint32_t> const& _info)
        : Ex{v[0]}
        , Ey{v[1]}
        , Ez{v[2]}
        , Bx{v[3]}
        , By{v[4]}
        , Bz{v[5]}
        , info{_info}
    {
        for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(states.size()); i++)
        {
            auto pos = i * 6;
            auto& em = states[i].electromag;
            Ex.send(em[0], _info[pos + 0]);
            Ey.send(em[1], _info[pos + 1]);
            Ez.send(em[2], _info[pos + 2]);
            Bx.send(em[3], _info[pos + 3]);
            By.send(em[4], _info[pos + 4]);
            Bz.send(em[5], _info[pos + 5]);
        }
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    Electromags() __device__
    {
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    EMInterop<Float> __device__ electromag(std::uint16_t i)
    {
        auto pos = i * 6;
        Electromags em;
        em.Ex   = this->Ex + this->info[pos + 0];
        em.Ey   = this->Ey + this->info[pos + 1];
        em.Ez   = this->Ez + this->info[pos + 2];
        em.Bx   = this->Bx + this->info[pos + 3];
        em.By   = this->By + this->info[pos + 4];
        em.Bz   = this->Bz + this->info[pos + 5];
        em.info = this->info + pos;
        return EMInterop<Float>{em};
    }

    container_t<Float> Ex, Ey, Ez, Bx, By, Bz;
    container_t<std::uint32_t> info;
};


template<typename SimOpts, bool GPU>
struct GridLayouts : kul::gpu::DeviceClass<GPU>
{
    static constexpr std::uint8_t dim  = SimOpts::dimension;
    static constexpr std::uint8_t dim2 = dim * 2;

    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = GridLayouts<SimOpts, true>;
    using Float = typename SimOpts::Float;

    using GridLayoutImpl = typename SimOpts::YeeLayout_t;
    using GridLayoutDAO  = PHARE::core::GridLayoutDAO<GridLayoutImpl, /*Refs=*/true>;
    using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl, GridLayoutDAO>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static auto make_shared(std::vector<PatchState<SimOpts>> const& states)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());
        return std::make_shared<GridLayouts<SimOpts, GPU>>(states, n_states);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    GridLayouts(std::vector<PatchState<SimOpts>> const& states, std::uint32_t n_states)
        : meshSize{n_states * dim}
        , origin{n_states * dim}
        , inverseMeshSize{n_states * dim}
        , nbrPhysicalCells{n_states * dim}
        , physicalStartIndices{n_states * dim * 2}
        , physicalEndIndices{n_states * dim * 2}
        , ghostEndIndices{n_states * dim * 2}
        , AMRBox{n_states * dim * 2}
    {
        for (std::uint32_t curr_part = 0, i = 0; i < n_states; i++)
        {
            GridLayout gridLayout{GridLayoutDAO{states[i].gridLayoutDAO}};
            this->add(gridLayout, i);
        }
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    void add(GridLayout const& gridLayout, std::uint32_t i)
    {
        std::uint32_t start = i * dim, start2 = i * dim2;

        meshSize.send(gridLayout.meshSize().data(), dim, start);
        origin.send(&gridLayout.origin()[0], dim, start);
        inverseMeshSize.send(gridLayout.inverseMeshSize().data(), dim, start);
        nbrPhysicalCells.send(gridLayout.nbrCells().data(), dim, start);
        physicalStartIndices.send(&gridLayout.physicalStartIndexTable()[0][0], dim2, start2);
        physicalEndIndices.send(&gridLayout.physicalEndIndexTable()[0][0], dim2, start2);
        ghostEndIndices.send(&gridLayout.ghostEndIndexTable()[0][0], dim2, start2);
        AMRBox.send(reinterpret_cast<std::int32_t const* const>(&gridLayout.AMRBox()), dim2,
                    start2);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(meshSize, origin, inverseMeshSize, nbrPhysicalCells,
                                            physicalStartIndices, physicalEndIndices,
                                            ghostEndIndices, AMRBox);
    }

    template<typename T>
    static auto __device__ array_cast(T* array)
    {
        return reinterpret_cast<std::array<T, dim>*>(array);
    }

    template<typename A, typename T>
    static auto __device__ array_cast(T* array)
    {
        return reinterpret_cast<A*>(array);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    GridLayout __device__ gridLayout(std::uint16_t i = 0)
    {
        constexpr std::uint32_t dim2 = dim * 2;
        std::uint32_t start = i * dim, start2 = i * dim2;
        return GridLayout{GridLayoutDAO{
            *array_cast(&meshSize[start]),
            *reinterpret_cast<core::Point<Float, dim>*>(array_cast(&origin[start])),
            *array_cast(&nbrPhysicalCells[start]), *array_cast(&inverseMeshSize[start]),
            *array_cast<std::array<std::array<std::uint32_t, dim>, 2>>(
                &physicalStartIndices[start2]),
            *array_cast<std::array<std::array<std::uint32_t, dim>, 2>>(&physicalEndIndices[start2]),
            *array_cast<std::array<std::array<std::uint32_t, dim>, 2>>(&ghostEndIndices[start2]),
            *reinterpret_cast<core::Box<std::int32_t, dim>*>(array_cast(&AMRBox[start2]))}};
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<Float> meshSize, origin, inverseMeshSize;
    container_t<std::uint32_t> nbrPhysicalCells;
    container_t<std::uint32_t> physicalStartIndices, physicalEndIndices, ghostEndIndices;
    container_t<std::int32_t> AMRBox; // <- upper+lower = dim * 2
};


template<typename SimOpts, bool GPU>
struct PatchStatePerParticle : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = PatchStatePerParticle<SimOpts, true>;
    using Float = typename SimOpts::Float;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    PatchStatePerParticle(std::uint32_t n_patches, std::uint32_t n_particles)
        : particles{1}
        , electromags{1}
        , gridLayouts{1}
        , particlePatchStateIdx{n_particles}
        , info{n_patches + 1}
    {
        info.send(&n_particles);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(particles, electromags, gridLayouts,
                                            particlePatchStateIdx, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ n_particles() const
    {
        return info[0];
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ operator[](std::size_t idx) const
    {
        return particlePatchStateIdx[idx];
    }

    container_t<EMContiguousParticles<Float, 1, true>> particles;
    container_t<Electromags<SimOpts, true>> electromags;
    container_t<GridLayouts<SimOpts, true>> gridLayouts;

    container_t<std::uint16_t> particlePatchStateIdx;
    container_t<std::uint32_t> info;
};

template<typename SimOpts>
struct ParticlePatchState
{
    static constexpr bool GPU         = false;
    static constexpr std::uint8_t dim = 1;
    using Float                       = typename SimOpts::Float;
    using GridLayout                  = typename GridLayouts<SimOpts, GPU>::GridLayout;
    using GridLayouts_                = GridLayouts<SimOpts, GPU>;
    using Electromags_                = Electromags<SimOpts, GPU>;
    using EMContiguousParticles_      = EMContiguousParticles<Float, dim, GPU>;
    using PatchStatePerParticle_      = PatchStatePerParticle<SimOpts, GPU>;

    ParticlePatchState(std::vector<PatchState<SimOpts>> const& states)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());

        for (auto const& data : states)
            for (auto const& particle_array : data.ions)
                n_particles += particle_array->size();

        particles = std::make_shared<EMContiguousParticles_>(n_particles);
        pspp      = std::make_shared<PatchStatePerParticle_>(n_states, n_particles);

        for (std::uint32_t parti = 0, i = 0; i < n_states; i++)
            for (auto const* particle_array : states[i].ions)
            {
                kul::gpu::fill_n(pspp->particlePatchStateIdx + parti, particle_array->size(), i);
                parti = particles->add(*particle_array, parti);
            }

        pspp->particles.send((*particles)());
        pspp->gridLayouts.send((*(gridLayouts = GridLayouts_::make_shared(states)))());
        pspp->electromags.send((*(electromags = Electromags_::make_shared(states)))());
    }

    auto operator()() { return (*pspp)(); }

    std::uint32_t n_particles = 0;
    std::shared_ptr<PatchStatePerParticle_> pspp;
    std::shared_ptr<EMContiguousParticles_> particles;
    std::shared_ptr<GridLayouts_> gridLayouts;
    std::shared_ptr<Electromags_> electromags;
};


static constexpr size_t X = 1024, Y = 1024, Z = 1;         // 40;
static constexpr size_t TPB_X = 16, TPB_Y = 16, TPB_Z = 1; // 4;
static constexpr size_t MAX_PARTICLES = X * Y * Z;

template<typename PHARE_TYPES>
__global__ void gpu_particles_in(PHARE::gpu::PatchStatePerParticle<PHARE_TYPES, true>* ppsp)
{
    auto i = kul::gpu::idx();
    if (i >= ppsp->n_particles())
        return;
    auto patchStateIDX         = (*ppsp)[i];
    auto gridLayout            = ppsp->gridLayouts->gridLayout(patchStateIDX);
    auto electromag            = ppsp->electromags->electromag(patchStateIDX);
    ppsp->particles->charge[i] = electromag.E.getComponent(PHARE::core::Component::X)(0);
}

template<typename Float, std::size_t dim = 1, std::size_t interp = 1, std::size_t nbRefineParts = 2>
void do_thing(std::string job_id)
{
    using PHARE_TYPES = PHARE::PHARE_Types<dim, interp, nbRefineParts, Float>;
    SimulatorTestParam<dim, interp, nbRefineParts, Float> sim{job_id};
    auto& hierarchy   = *sim.hierarchy;
    auto& hybridModel = *sim.getHybridModel();
    auto& state       = hybridModel.state;
    auto topLvl       = hierarchy.getNumberOfLevels() - 1;

    std::vector<PHARE::gpu::PatchState<PHARE_TYPES>> states;
    PHARE::amr::visitHierarchy<typename PHARE_TYPES::GridLayout_t>(
        hierarchy, *hybridModel.resourcesManager,
        [&](auto& gridLayout, std::string, size_t) { states.emplace_back(gridLayout, state); },
        topLvl, topLvl + 1, hybridModel);

    PHARE::gpu::ParticlePatchState<PHARE_TYPES> packer{states};
    kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(gpu_particles_in<PHARE_TYPES>, packer());

    KLOG(NON) << "MAX_PARTICLES: " << MAX_PARTICLES;
    KLOG(NON) << "GPU PARTICLES: " << packer.n_particles;

    for (auto const& state : states)
        KLOG(NON) << state.electromag[0][1];

    auto charge = packer.particles->charge();
    KLOG(NON) << charge[0];
    KLOG(NON) << charge.back();
}

} // namespace PHARE::gpu

#endif /*PHARE_GPU_HPP*/
