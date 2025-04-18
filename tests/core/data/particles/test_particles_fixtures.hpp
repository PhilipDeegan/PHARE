#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP

#include "test_particles.hpp"

#include <string>
#include <cassert>

namespace PHARE::core
{

template<typename ParticleArray_t>
struct UsableParticlesPopulation
{
    template<typename GridLayout_t>
    UsableParticlesPopulation(std::string const& _name, GridLayout_t const& layout)
        : name{_name}
        , domain_particles{make_particles<ParticleArray_t>(layout)}
        , patch_ghost_particles{make_particles<ParticleArray_t>(layout)}
        , level_ghost_particles{make_particles<ParticleArray_t>(layout)}
    {
    }

    auto& pack() { return particles_pack; }
    auto& pack() const { return particles_pack; }

    std::string name;
    ParticleArray_t domain_particles;
    ParticleArray_t patch_ghost_particles;
    ParticleArray_t level_ghost_particles;
    core::ParticlesPack<ParticleArray_t> particles_pack{name,
                                                        &domain_particles,
                                                        &patch_ghost_particles,
                                                        &level_ghost_particles,
                                                        /*levelGhostParticlesOld=*/nullptr,
                                                        /*levelGhostParticlesNew=*/nullptr};
};


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
