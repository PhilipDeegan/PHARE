#ifndef PHARE_GPU_HPP
#define PHARE_GPU_HPP

#include "kul/gpu/rocm.hpp"

#include "tests/simulator/per_test.h"
using PHARE_TYPES = PHARE::PHARE_Types<dim, interp, nbRefineParts>;



namespace PHARE::gpu
{
template<size_t dim, bool GPU>
struct EMContiguousParticles : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = EMContiguousParticles<dim, true>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    EMContiguousParticles(size_t nbr) // bytes per particle
        : leaving{nbr}                // 1
        , iCell{nbr * dim}            // 4 * dim
        , delta{nbr * dim}            // 4 * dim
        , weight{nbr}                 // 8
        , charge{nbr}                 // 8
        , v{nbr * 3}                  // 24
        , E{nbr * 3}                  // 24
        , B{nbr * 3}                  // 24
    {                                 // 1d  97 - 1e6 ~100MB / 1e7 ~1GB
    }                                 // 2d 105
                                      // 3d 113

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) add(core::ParticleArray<dim> const& array, size_t i)
    {
        PHARE::core::EMContiguousParticles<dim> particles{array};
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
    decltype(auto) operator()()
    {
        return Super::template alloc<gpu_t>(leaving, iCell, delta, weight, charge, v, E, B);
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<bool> leaving;
    container_t<int32_t> iCell;
    container_t<float> delta;
    container_t<double> weight, charge, v, E, B;
};

template<bool GPU>
struct Electromag : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = Electromag<true>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    Electromag(std::vector<kul::Pointers<double>> const& v)
        : Ex{v[0]}
        , Ey{v[1]}
        , Ez{v[2]}
        , Bx{v[3]}
        , By{v[4]}
        , Bz{v[5]}
    {
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz);
    }

    typename Super::template container_t<double> Ex, Ey, Ez, Bx, By, Bz;
};

template<bool GPU>
struct GridLayouts : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = GridLayouts<true>;

    using GridLayoutImpl      = PHARE_TYPES::YeeLayout_t;
    static constexpr auto dim = GridLayoutImpl::dimension;
    using GridLayoutDAO       = PHARE::core::GridLayoutDAO<dim, /*Refs=*/true>;
    using GridLayout          = PHARE::core::GridLayout<GridLayoutImpl, GridLayoutDAO>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    GridLayouts(size_t size)
        : meshSize{size * dim}
        , origin{size * dim}
        , inverseMeshSize{size * dim}
        , nbrPhysicalCells{size * dim}
        , physicalStartIndices{size * dim * 2}
        , physicalEndIndices{size * dim * 2}
        , ghostEndIndices{size * dim * 2}
        , AMRBox{size * dim * 2}
    {
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    void add(GridLayout const& gridLayout, size_t i)
    {
        constexpr size_t dim2 = dim * 2;
        size_t start = i * dim, start2 = i * dim2;

        meshSize.send(gridLayout.meshSize().data(), dim, start);
        origin.send(&gridLayout.origin()[0], dim, start);
        inverseMeshSize.send(gridLayout.inverseMeshSize().data(), dim, start);
        nbrPhysicalCells.send(gridLayout.nbrCells().data(), dim, start);
        physicalStartIndices.send(&gridLayout.physicalStartIndexTable()[0][0], dim2, start2);
        physicalEndIndices.send(&gridLayout.physicalEndIndexTable()[0][0], dim2, start2);
        ghostEndIndices.send(&gridLayout.ghostEndIndexTable()[0][0], dim2, start2);
        AMRBox.send(reinterpret_cast<int32_t const* const>(&gridLayout.AMRBox()), dim2, start2);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) operator()()
    {
        return Super::template alloc<gpu_t>(meshSize, origin, inverseMeshSize, nbrPhysicalCells,
                                            physicalStartIndices, physicalEndIndices,
                                            ghostEndIndices, AMRBox);
    }

    template<typename T>
    static decltype(auto) array_cast(T* array)
    {
        return reinterpret_cast<std::array<T, dim>*>(array);
    }

    template<typename A, typename T>
    static decltype(auto) array_cast(T* array)
    {
        return reinterpret_cast<A*>(array);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto gridLayout(uint16_t i = 0)
    {
        constexpr size_t dim2 = dim * 2;
        size_t start = i * dim, start2 = i * dim2;
        return GridLayout{GridLayoutDAO{
            *array_cast(&meshSize[start]),
            *reinterpret_cast<core::Point<double, dim>*>(array_cast(&origin[start])),
            *array_cast(&nbrPhysicalCells[start]), *array_cast(&inverseMeshSize[start]),
            *array_cast<std::array<std::array<uint32_t, dim>, 2>>(&physicalStartIndices[start2]),
            *array_cast<std::array<std::array<uint32_t, dim>, 2>>(&physicalEndIndices[start2]),
            *array_cast<std::array<std::array<uint32_t, dim>, 2>>(&ghostEndIndices[start2]),
            *reinterpret_cast<core::Box<int32_t, dim>*>(array_cast(&AMRBox[start2]))}};
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<double> meshSize, origin, inverseMeshSize;
    container_t<uint32_t> nbrPhysicalCells;
    container_t<uint32_t> physicalStartIndices, physicalEndIndices, ghostEndIndices;
    container_t<int32_t> AMRBox; // <- upper+lower = dim * 2
};



struct PatchState
{
    using GridLayout = PHARE_TYPES::GridLayout_t;

    template<typename State>
    PatchState(GridLayout const& gridLayout, State const& state)
        : gridLayoutDAO{gridLayout}
    {
        for (auto const& pop : state.ions)
            ions.emplace_back(&pop.domainParticles());

        auto vecF = [&](auto& EBxyz) {
            for (size_t i = 0; i < 3; i++)
                electromag.emplace_back(EBxyz[i].data(), EBxyz[i].size());
        };

        vecF(state.electromag.E);
        vecF(state.electromag.B);
    }

    GridLayout::Super gridLayoutDAO;
    std::vector<kul::Pointers<double>> electromag;
    std::vector<core::ParticleArray<dim>*> ions;
};


template<bool GPU>
struct PatchStatePerParticle : kul::gpu::DeviceClass<GPU>
{
    using gpu_t = PatchStatePerParticle<true>;
    using Super = kul::gpu::DeviceClass<GPU>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    PatchStatePerParticle(size_t n_patches, size_t n_particles)
        : particles{1}
        , electromags{n_patches}
        , gridLayouts{1}
        , patchStatePerParticle{n_particles}
        , info{n_patches + 1}
    {
        info.send(&n_particles);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) operator()()
    {
        return Super::template alloc<gpu_t>(particles, electromags, gridLayouts,
                                            patchStatePerParticle, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    size_t n_particles() const
    {
        return info[0];
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<EMContiguousParticles<1, true>> particles;
    container_t<Electromag<true>> electromags;
    container_t<GridLayouts<true>> gridLayouts;
    container_t<uint16_t> patchStatePerParticle;
    container_t<size_t> info;
};

struct ParticlePatchState
{
    static constexpr bool GPU = false;
    using GridLayout          = GridLayouts<false>::GridLayout;

    ParticlePatchState(std::vector<PatchState> const& v)
    {
        for (auto const& data : v)
            for (auto const& particle_array : data.ions)
                n_particles += particle_array->size();

        this->particles   = std::make_shared<EMContiguousParticles<1, GPU>>(n_particles);
        this->gridLayouts = std::make_shared<GridLayouts<false>>(v.size());

        statePack  = std::make_shared<PatchStatePerParticle<GPU>>(v.size(), n_particles);
        auto& pack = *statePack;
        for (size_t curr_part = 0, i = 0; i < v.size(); i++)
        {
            for (auto const* particle_array : v[i].ions)
            {
                kul::gpu::fill_n(pack.patchStatePerParticle + curr_part, particle_array->size(), i);
                curr_part = this->particles->add(*particle_array, curr_part);
            }
            gpu_electromags.emplace_back((*this->electromags.emplace_back(
                std::make_shared<Electromag<GPU>>(v[i].electromag)))());
            GridLayout gridLayout{GridLayout::Super{v[i].gridLayoutDAO}};
            this->gridLayouts->add(gridLayout, i);
        }
        pack.particles.send((*particles)());
        pack.gridLayouts.send((*gridLayouts)());
        pack.electromags.send(gpu_electromags);
    }

    decltype(auto) operator()() { return (*statePack)(); }

    size_t n_particles = 0;
    std::shared_ptr<PatchStatePerParticle<GPU>> statePack;
    std::shared_ptr<EMContiguousParticles<1, GPU>> particles;
    std::shared_ptr<GridLayouts<GPU>> gridLayouts;
    std::vector<std::shared_ptr<Electromag<GPU>>> electromags;
    std::vector<Electromag<true>*> gpu_electromags;
};


__global__ void particles_in(PatchStatePerParticle<true>* ppsp)
{
    auto i = kul::gpu::hip::idx();
    if (i > ppsp->n_particles())
        return;
    auto patchStateIDX         = ppsp->patchStatePerParticle[i];
    auto gridLayout            = ppsp->gridLayouts->gridLayout(patchStateIDX);
    auto electromag            = ppsp->electromags[patchStateIDX];
    ppsp->particles->charge[i] = gridLayout.ghostEndIndexTable()[0][0];
}

} // namespace PHARE::gpu


#endif /*PHARE_GPU_HPP*/