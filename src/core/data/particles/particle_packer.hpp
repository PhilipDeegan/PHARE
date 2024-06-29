#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP

#include <vector>
#include <cstddef>

#include "core/def.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"


namespace PHARE::core
{
// PGI compiler (nvc++ 21.3-0) doesn't like static initializations of arrays,
//   would result in empty strings
inline std::array<std::string, 5> packer_keys()
{
    // The order of this array must match the tuple order of ParticlePacker::get(particle)
    return {"weight", "charge", "iCell", "delta", "v"};
}

template<typename ParticleArray_>
class ParticlePacker
{
    auto constexpr static dim = ParticleArray_::dimension;
    using ParticleArray_t     = ParticleArray_;
    using Particle_t          = core::Particle<dim>;

public:
    static constexpr std::size_t n_keys = 5;

    ParticlePacker(ParticleArray_t const& particles)
        : particles_{particles}
    {
    }


    NO_DISCARD static auto get(Particle_t const& particle)
    {
        return std::forward_as_tuple(particle.weight_, particle.charge_, particle.iCell_,
                                     particle.delta_, particle.v_);
    }

    static auto empty()
    {
        Particle<dim> particle{};
        return get(particle);
    }


    // sometimes we use this to infer the size of an ParticleArray
    // could be "charge" either
    NO_DISCARD static auto arbitrarySingleValueKey() { return "weight"; }

    NO_DISCARD static auto keys() { return packer_keys(); }

    NO_DISCARD auto get(std::size_t i) const { return get(particles_[i]); }
    NO_DISCARD bool hasNext() const { return it_ < particles_.size(); }
    NO_DISCARD auto next() { return get(it_++); }

    template<typename SoAParticles_t>
    void pack(SoAParticles_t& copy)
    {
        if constexpr (ParticleArray_t::layout_mode == LayoutMode::SoA)
            copy = particles_;
        else
        {
            std::size_t idx = 0;
            while (this->hasNext())
            {
                auto next        = this->next();
                copy.weight(idx) = std::get<0>(next);
                copy.charge(idx) = std::get<1>(next);
                copy.iCell(idx)  = std::get<2>(next);
                copy.delta(idx)  = std::get<3>(next);
                copy.v(idx)      = std::get<4>(next);
                ++idx;
            }
        }
    }

private:
    ParticleArray_t const& particles_;
    std::size_t it_ = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
