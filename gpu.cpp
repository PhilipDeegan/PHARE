
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

static constexpr size_t X = 1024, Y = 1024, Z = 40;
static constexpr size_t TPB_X = 16, TPB_Y = 16, TPB_Z = 4;
static constexpr size_t MAX_PARTICLES = X * Y * Z;
static constexpr size_t dim = 1, interp = 1, nbRefineParts = 2;

#include "gpu.hpp"

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    TestSimulator<dim, interp> sim;
    auto& hierarchy   = *sim.hierarchy;
    auto& hybridModel = *sim.getHybridModel();
    auto& state       = hybridModel.state;
    auto topLvl       = hierarchy.getNumberOfLevels() - 1;

    std::vector<PHARE::gpu::PatchState> states;

    PHARE::amr::visitHierarchy<PHARE_TYPES::GridLayout_t>(
        hierarchy, *hybridModel.resourcesManager,
        [&](auto& gridLayout, std::string, size_t) { states.emplace_back(gridLayout, state); },
        topLvl, topLvl + 1, hybridModel);

    PHARE::gpu::ParticlePatchState packer{states};

    kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(PHARE::gpu::particles_in, packer());

    KLOG(NON) << "MAX_PARTICLES: " << MAX_PARTICLES;
    KLOG(NON) << "GPU PARTICLES: " << packer.n_particles;
    KLOG(NON) << packer.particles->charge()[packer.n_particles - 1];

    return 0;
}


/*

__global__ void xyz(uint32_t* in)
{
    auto i = kul::gpu::hip::idx();
    in[i]  = i;
}
int main(int argc, char** argv)
{
    std::vector<uint32_t> data(X * Y * Z, 123);
    kul::gpu::DeviceMem<uint32_t> gata{data};

    kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(xyz, gata.p);

    auto host = gata();

    KLOG(NON) << host[10];
    KLOG(NON) << host[100];
    KLOG(NON) << host[1000];
    KLOG(NON) << host[10000];
    KLOG(NON) << host[100000];
    KLOG(NON) << host[1000000];
    KLOG(NON) << host[2000000];
    KLOG(NON) << host.back();
    KLOG(NON) << host.size();

    return 0;
}*/