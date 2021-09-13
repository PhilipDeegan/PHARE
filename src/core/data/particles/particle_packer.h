#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_H
#define PHARE_CORE_DATA_PARTICLE_PACKER_H

#include <vector>
#include <cstddef>

#include "particle_array.h"
#if defined(HAVE_UMPIRE)
#include "llnl/particle_array.h"
#endif

namespace PHARE::core
{
// PGI compiler (nvc++ 21.3-0) doesn't like static initializations of arrays,
//   would result in empty strings
inline std::array<std::string, 5> packer_keys()
{
    // The order of this array must match the tuple order of ParticlePacker::get(particle)
    return {"weight", "charge", "iCell", "delta", "v"};
}

template<typename Particle>
class ParticlePacker
{
    static constexpr auto dim = Particle::dimension;

public:
    ParticlePacker(ParticleArray<Particle> const& particles)
        : particles_{particles}
    {
    }

#if defined(HAVE_UMPIRE)
    auto to_host_particle_array(llnl::ParticleArray<Particle> const& particles) {}

    ParticlePacker(llnl::ParticleArray<Particle> const& particles)
        : particles_{particles()}
    {
    }
#endif

    static auto get(Particle const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static auto empty()
    {
        Particle particle;
        return get(particle);
    }


    auto get(std::size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

    void pack(ContiguousParticles<dim>& copy)
    {
        auto copyTo = [](auto& a, auto& idx, auto size, auto& v) {
            std::copy(a.begin(), a.begin() + size, v.begin() + (idx * size));
        };
        std::size_t idx = 0;
        while (this->hasNext())
        {
            auto next        = this->next();
            copy.weight[idx] = std::get<0>(next);
            copy.charge[idx] = std::get<1>(next);
            copyTo(std::get<2>(next), idx, dim, copy.iCell);
            copyTo(std::get<3>(next), idx, dim, copy.delta);
            copyTo(std::get<4>(next), idx, 3, copy.v);
            idx++;
        }
    }

private:
    ParticleArray<Particle> const& particles_;
    std::size_t it_ = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
