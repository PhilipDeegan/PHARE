#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_H
#define PHARE_CORE_DATA_PARTICLE_PACKER_H

#include <cstddef>

#include "core/data/particles/particle.h"
#include "core/data/particles/particle_array.h"
#include "core/data/particles/particle_array_soa.h"

namespace PHARE::core
{
template<std::size_t dim>
class ParticlePacker
{
public:
    ParticlePacker(ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    static auto get(Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight_, particle.charge_, particle.iCell_,
                                     particle.delta_, particle.v_);
    }

    static auto empty()
    {
        Particle<dim> particle;
        return get(particle);
    }

    static auto& keys() { return keys_; }

    auto get(std::size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

    template<typename ParticleArray_SOA_t>
    void pack(ParticleArray_SOA_t& copy)
    {
        auto copyTo = [](auto& a, auto& idx, auto size, auto& v) {
            std::copy(a.begin(), a.begin() + size, &v[idx]);
        };
        std::size_t idx = 0;
        while (this->hasNext())
        {
            auto next        = this->next();
            copy.weight[idx] = std::get<0>(next);
            copy.charge[idx] = std::get<1>(next);
            copy.iCell_[idx] = std::get<2>(next);
            copy.delta[idx]  = std::get<3>(next);
            copy.v[idx]      = std::get<4>(next);
            idx++;
        }
    }

private:
    ParticleArray<dim> const& particles_;
    std::size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};



template<std::size_t dim>
ParticleArray_SOA<dim>::ParticleArray_SOA(ParticleArray<dim> const& particles)
    : ParticleArray_SOA{particles.size()}
{
    ParticlePacker<dim>{particles}.pack(*this);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
