#ifndef PHARE_GPU_HPP
#define PHARE_GPU_HPP

#include "kul/gpu/rocm.hpp"

#include "tests/simulator/per_test.h"
using PHARE_TYPES = PHARE::PHARE_Types<dim, interp, nbRefineParts>;



namespace PHARE::gpu
{
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
            for (uint8_t i = 0; i < 3; i++)
                electromag.emplace_back(EBxyz[i].data(), EBxyz[i].size());
        };

        vecF(state.electromag.E);
        vecF(state.electromag.B);
    }

    GridLayout::Super gridLayoutDAO;
    std::vector<kul::Pointers<double, uint32_t>> electromag;
    std::vector<core::ParticleArray<dim>*> ions;
};


template<uint8_t dim, bool GPU>
struct EMContiguousParticles : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = EMContiguousParticles<dim, true>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    EMContiguousParticles(uint32_t nbr) // bytes per particle
        : leaving{nbr}                  // 1
        , iCell{nbr * dim}              // 4 * dim
        , delta{nbr * dim}              // 4 * dim
        , weight{nbr}                   // 8
        , charge{nbr}                   // 8
        , v{nbr * 3}                    // 24
        , E{nbr * 3}                    // 24
        , B{nbr * 3}                    // 24
    {                                   // 1d  97 - 1e6 ~100MB / 1e7 ~1GB
    }                                   // 2d 105
                                        // 3d 113

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) add(core::ParticleArray<dim> const& array, uint32_t i)
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


struct FieldInterop
{
    using Pointers = kul::Pointers<double, uint32_t>;

    decltype(auto) operator()(uint32_t i) { return ptrs[i]; }
    // decltype(auto) operator()(uint32_t, uint32_t) { return 2; }
    // decltype(auto) operator()(uint32_t, uint32_t, uint32_t) { return 3; }

    Pointers ptrs;
};

struct VecFieldInterop
{
    using Pointers = FieldInterop::Pointers;

    template<typename Electromags>
    static VecFieldInterop E(Electromags const& em)
    {
        return VecFieldInterop{
            {{{em.Ex, em.info[0]}}, {{em.Ey, em.info[1]}}, {{em.Ez, em.info[2]}}}};
    }
    template<typename Electromags>
    static VecFieldInterop B(Electromags const& em)
    {
        return VecFieldInterop{
            {{{em.Bx, em.info[3]}}, {{em.By, em.info[4]}}, {{em.Bz, em.info[5]}}}};
    }


    FieldInterop& getComponent(PHARE::core::Component XYZ)
    {
        if (XYZ == PHARE::core::Component::X)
            return ptrs[0];
        else if (XYZ == PHARE::core::Component::Y)
            return ptrs[1];

        return ptrs[2];
    }

    FieldInterop ptrs[3];
};

struct EMInterop
{
    template<typename Electromags>
    EMInterop(Electromags&& _em)
        : E{VecFieldInterop::E(_em)}
        , B{VecFieldInterop::B(_em)}
    {
    }

    VecFieldInterop E, B;
};


template<bool GPU>
struct Electromags : kul::gpu::DeviceClass<GPU>
{
    using Super    = kul::gpu::DeviceClass<GPU>;
    using Pointers = kul::Pointers<double, uint32_t>;
    using gpu_t    = Electromags<true>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static decltype(auto) make_shared(std::vector<PatchState> const& states)
    {
        uint32_t n_states = static_cast<uint32_t>(states.size());

        std::vector<uint32_t> emXYZ(6, 0), emInfo(n_states * 12);
        for (size_t j = 0; j < n_states; j++)
            for (size_t i = 0; i < states[j].electromag.size(); i++)
            {
                auto pos        = (j * 12) + i;
                emInfo[pos]     = states[j].electromag[i].size();
                emInfo[pos + 6] = emXYZ[i];
                emXYZ[i] += emInfo[pos];
            }

        return std::make_shared<Electromags<GPU>>(states, emXYZ, emInfo);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    Electromags(std::vector<PatchState> const& states, std::vector<uint32_t> const& v,
                std::vector<uint32_t> const& _info)
        : Ex{v[0]}
        , Ey{v[1]}
        , Ez{v[2]}
        , Bx{v[3]}
        , By{v[4]}
        , Bz{v[5]}
        , info{_info}
    {
        for (uint32_t i = 0; i < static_cast<uint32_t>(states.size()); i++)
        {
            auto pos = i * 12;
            auto& em = states[i].electromag;
            Ex.send(em[0], _info[pos + 0 + 6]);
            Ey.send(em[1], _info[pos + 1 + 6]);
            Ez.send(em[2], _info[pos + 2 + 6]);
            Bx.send(em[3], _info[pos + 3 + 6]);
            By.send(em[4], _info[pos + 4 + 6]);
            Bz.send(em[5], _info[pos + 5 + 6]);
        }
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    decltype(auto) operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    Electromags()
    {
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    decltype(auto) electromag(uint16_t i)
    {
        auto pos = i * 12;
        Electromags em;
        em.Ex   = this->Ex + this->info[pos + 0 + 6];
        em.Ey   = this->Ey + this->info[pos + 1 + 6];
        em.Ez   = this->Ez + this->info[pos + 2 + 6];
        em.Bx   = this->Bx + this->info[pos + 3 + 6];
        em.By   = this->By + this->info[pos + 4 + 6];
        em.Bz   = this->Bz + this->info[pos + 5 + 6];
        em.info = this->info + pos;
        return EMInterop{em};
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<double> Ex, Ey, Ez, Bx, By, Bz;
    container_t<uint32_t> info;
};


template<bool GPU>
struct GridLayouts : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = GridLayouts<true>;

    using GridLayoutImpl          = PHARE_TYPES::YeeLayout_t;
    static constexpr uint8_t dim  = GridLayoutImpl::dimension;
    static constexpr uint8_t dim2 = dim * 2;
    using GridLayoutDAO           = PHARE::core::GridLayoutDAO<dim, /*Refs=*/true>;
    using GridLayout              = PHARE::core::GridLayout<GridLayoutImpl, GridLayoutDAO>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static decltype(auto) make_shared(std::vector<PatchState> const& states)
    {
        uint32_t n_states = static_cast<uint32_t>(states.size());
        return std::make_shared<GridLayouts<GPU>>(states, n_states);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    GridLayouts(std::vector<PatchState> const& states, uint32_t n_states)
        : meshSize{n_states * dim}
        , origin{n_states * dim}
        , inverseMeshSize{n_states * dim}
        , nbrPhysicalCells{n_states * dim}
        , physicalStartIndices{n_states * dim * 2}
        , physicalEndIndices{n_states * dim * 2}
        , ghostEndIndices{n_states * dim * 2}
        , AMRBox{n_states * dim * 2}
    {
        for (uint32_t curr_part = 0, i = 0; i < n_states; i++)
        {
            GridLayout gridLayout{GridLayoutDAO{states[i].gridLayoutDAO}};
            this->add(gridLayout, i);
        }
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    void add(GridLayout const& gridLayout, uint32_t i)
    {
        uint32_t start = i * dim, start2 = i * dim2;

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
    decltype(auto) gridLayout(uint16_t i = 0)
    {
        constexpr uint32_t dim2 = dim * 2;
        uint32_t start = i * dim, start2 = i * dim2;
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


template<bool GPU>
struct PatchStatePerParticle : kul::gpu::DeviceClass<GPU>
{
    using gpu_t = PatchStatePerParticle<true>;
    using Super = kul::gpu::DeviceClass<GPU>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    PatchStatePerParticle(uint32_t n_patches, uint32_t n_particles)
        : particles{1}
        , electromags{1}
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
    uint32_t n_particles() const
    {
        return info[0];
    }

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    container_t<EMContiguousParticles<1, true>> particles;
    container_t<Electromags<true>> electromags;
    container_t<GridLayouts<true>> gridLayouts;

    container_t<uint16_t> patchStatePerParticle;
    container_t<uint32_t> info;
};

struct ParticlePatchState
{
    static constexpr bool GPU = false;
    using GridLayout          = GridLayouts<false>::GridLayout;

    ParticlePatchState(std::vector<PatchState> const& states)
    {
        uint32_t n_states = static_cast<uint32_t>(states.size());

        for (auto const& data : states)
            for (auto const& particle_array : data.ions)
                n_particles += particle_array->size();

        this->particles = std::make_shared<EMContiguousParticles<1, GPU>>(n_particles);
        this->statePack = std::make_shared<PatchStatePerParticle<GPU>>(n_states, n_particles);
        auto& pack      = *this->statePack;

        for (uint32_t curr_part = 0, i = 0; i < n_states; i++)
            for (auto const* particle_array : states[i].ions)
            {
                kul::gpu::fill_n(pack.patchStatePerParticle + curr_part, particle_array->size(), i);
                curr_part = this->particles->add(*particle_array, curr_part);
            }

        pack.particles.send((*particles)());
        pack.gridLayouts.send((*(this->gridLayouts = GridLayouts<GPU>::make_shared(states)))());
        pack.electromags.send((*(this->electromags = Electromags<GPU>::make_shared(states)))());
    }

    decltype(auto) operator()() { return (*statePack)(); }

    uint32_t n_particles = 0;
    std::shared_ptr<PatchStatePerParticle<GPU>> statePack;
    std::shared_ptr<EMContiguousParticles<1, GPU>> particles;
    std::shared_ptr<GridLayouts<GPU>> gridLayouts;
    std::shared_ptr<Electromags<GPU>> electromags;
};


__global__ void particles_in(PatchStatePerParticle<true>* ppsp)
{
    auto i = kul::gpu::hip::idx();
    if (i > ppsp->n_particles())
        return;
    auto patchStateIDX = ppsp->patchStatePerParticle[i];
    auto gridLayout    = ppsp->gridLayouts->gridLayout(patchStateIDX);
    auto electromag    = ppsp->electromags->electromag(patchStateIDX);

    ppsp->particles->charge[i] = electromag.E.getComponent(PHARE::core::Component::X)(1);

    // ppsp->particles->charge[i] = electromag.Ex[2];

    // ppsp->particles->charge[i] = gridLayout.ghostEndIndexTable()[0][0];

    // if (patchStateIDX == 0)
    // ppsp->particles->charge[i] = emFields.E.getComponent(PHARE::core::Component::X)(1);
    // electromag.info[0]; // electromag.Ex[1];
    // else
    //     ppsp->particles->charge[i] = patchStateIDX; // gridLayout.ghostEndIndexTable()[0][0];
}

} // namespace PHARE::gpu


#endif /*PHARE_GPU_HPP*/