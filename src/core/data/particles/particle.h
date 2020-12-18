#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include <array>
#include <random>
#include <type_traits>

#include "core/utilities/point/point.h"
#include "core/utilities/span.h"
#include "core/utilities/types.h"


namespace PHARE::core
{
template<typename T = float>
struct ParticleDeltaDistribution
{
    constexpr static T epsilon = std::numeric_limits<T>::epsilon();
    constexpr static T one     = 1;

    template<typename Generator>
    float operator()(Generator& generator)
    {
        return dist(generator);
    }
    std::uniform_real_distribution<T> dist{0, one - epsilon};
};


template<typename Particle>
auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}



template<typename Float, size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");

    using float_type              = Float;
    static const size_t dimension = dim;

    Float weight;
    Float charge;

    std::array<int, dim> iCell   = ConstArray<int, dim>();
    std::array<float, dim> delta = ConstArray<float, dim>();
    std::array<Float, 3> v       = ConstArray<Float, 3>();

    Float Ex = 0, Ey = 0, Ez = 0;
    Float Bx = 0, By = 0, Bz = 0;
};




template<typename Float, std::size_t dim>
struct ParticleView
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");

    using float_type                       = Float;
    static constexpr std::size_t dimension = dim;

    Float& weight;
    Float& charge;
    std::array<int, dim>& iCell;
    std::array<float, dim>& delta;
    std::array<Float, 3>& v;
};



template<typename Float, std::size_t dim, typename T>
inline constexpr auto is_phare_particle_type
    = std::is_same_v<Particle<Float, dim>, T> or std::is_same_v<ParticleView<Float, dim>, T>;


template<typename Float, std::size_t dim, template<typename, std::size_t> typename ParticleA,
         template<typename, std::size_t> typename ParticleB>
typename std::enable_if_t<                                    //
    is_phare_particle_type<Float, dim, ParticleA<Float, dim>> //
        and is_phare_particle_type<Float, dim, ParticleB<Float, dim>>,
    bool>
operator==(ParticleA<Float, dim> const& particleA, ParticleB<Float, dim> const& particleB)
{
    return particleA.weight == particleB.weight and //
           particleA.charge == particleB.charge and //
           particleA.iCell == particleB.iCell and   //
           particleA.delta == particleB.delta and   //
           particleA.v == particleB.v;
}

} // namespace PHARE::core


namespace std
{
template<typename Float, size_t dim, template<typename, std::size_t> typename Particle_t>
typename std::enable_if_t<PHARE::core::is_phare_particle_type<Float, dim, Particle_t<Float, dim>>,
                          PHARE::core::Particle<Float, dim>>
copy(Particle_t<Float, dim> const& from)
{
    return {from.weight, from.charge, from.iCell, from.delta, from.v};
}

} // namespace std


#endif
