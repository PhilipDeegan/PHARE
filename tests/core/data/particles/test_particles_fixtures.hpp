#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP

#include <cassert>
#include <functional>

#include "phare_core.hpp"
#include "core/utilities/types.hpp"


namespace PHARE::core
{

template<typename ParticleArray_t>
struct UsableParticlesPopulation
{
    template<typename... Args>
    UsableParticlesPopulation(Args&&... args)
        : domain_particles{args...}
    {
    }

    auto pack() { return &particles_pack; }

    ParticleArray_t domain_particles;
    ParticleArray_t patch_ghost_particles = domain_particles;
    ParticleArray_t level_ghost_particles = domain_particles;
    core::ParticlesPack<ParticleArray_t> particles_pack{&domain_particles, //
                                                        &patch_ghost_particles,
                                                        &level_ghost_particles,
                                                        /*levelGhostParticlesOld=*/nullptr,
                                                        /*levelGhostParticlesNew=*/nullptr};
};



template<std::size_t dim>
PHARE::core::Particle<dim> particle(std::array<int, dim> const& icell)
{
    return {/*.weight = */ .001,
            /*.charge = */ .001,
            /*.iCell  = */ icell,
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.505, .505, .505}}};
}

template<std::size_t dim>
PHARE::core::Particle<dim> particle(int const icell = 15)
{
    return particle(PHARE::core::ConstArray<int, dim>(icell));
}

template<typename Particles>
void disperse(Particles& particles, std::optional<int> seed = std::nullopt)
{
    auto gen = [&]() {
        if (seed.has_value())
            return std::mt19937_64(*seed);
        std::random_device rd;
        std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd()};
        return std::mt19937_64(seed_seq);
    }();

    std::shuffle(particles.begin(), particles.end(), gen);

    // reset cellmap
    particles.empty_map();
    particles.map_particles();
}

template<typename Particles, std::size_t ghost_cells = 0, typename GridLayout_t>
auto make_particles(GridLayout_t const& layout)
{
    return Particles{grow(layout.AMRBox(), ghost_cells + 1)};
}


template<typename Particles, typename Box>
auto add_particles_in(Particles& particles, Box const& box, std::size_t const ppc)
{
    for (auto const& iterator : box)
        for (std::size_t i = 0; i < ppc; ++i)
            particles.push_back(particle(iterator.template toArray()));

    return particles;
}



} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
