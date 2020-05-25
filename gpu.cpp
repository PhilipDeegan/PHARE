
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

static constexpr size_t X = 1024, Y = 1024, Z = 40;
static constexpr size_t TPB_X = 16, TPB_Y = 16, TPB_Z = 4;
static constexpr size_t MAX_PARTICLES = X * Y * Z;
static constexpr size_t dim = 1, interp = 1, nbRefineParts = 2;

#include "gpu.hpp"

template<typename Float>
void do_thing()
{
    using PHARE_TYPES = PHARE::PHARE_Types<dim, interp, nbRefineParts, Float>;
    using SimOpts     = PHARE_TYPES;
    TestSimulator<dim, interp, nbRefineParts, Float> sim;
    auto& hierarchy   = *sim.hierarchy;
    auto& hybridModel = *sim.getHybridModel();
    auto& state       = hybridModel.state;
    auto topLvl       = hierarchy.getNumberOfLevels() - 1;

    std::vector<PHARE::gpu::PatchState<SimOpts>> states;
    PHARE::amr::visitHierarchy<typename PHARE_TYPES::GridLayout_t>(
        hierarchy, *hybridModel.resourcesManager,
        [&](auto& gridLayout, std::string, size_t) { states.emplace_back(gridLayout, state); },
        topLvl, topLvl + 1, hybridModel);

    for (auto const& state : states)
        KLOG(OTH) << state.electromag[0].p[1];

    PHARE::gpu::ParticlePatchState<SimOpts> packer{states};

    kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(PHARE::gpu::particles_in<SimOpts>, packer());

    KLOG(NON) << "MAX_PARTICLES: " << MAX_PARTICLES;
    KLOG(NON) << "GPU PARTICLES: " << packer.n_particles;
    auto charge = packer.particles->charge();

    KLOG(NON) << charge[0];
    KLOG(NON) << charge.back();

    for (auto const& state : states)
        KLOG(OTH) << state.electromag[0].p[1];
}

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    // if (argc == 1)
    //     do_thing<double>();
    // else
    do_thing<float>();
    return 0;
}
