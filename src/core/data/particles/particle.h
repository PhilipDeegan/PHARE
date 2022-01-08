#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include <array>
#include <random>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <iostream>

#include "core/utilities/point/point.h"
#include "core/utilities/span.h"
#include "core/utilities/types.h"


namespace PHARE::core
{
template<typename T = float>
struct ParticleDeltaDistribution
{
    template<typename Generator>
    T operator()(Generator& generator)
    {
        return dist(generator);
    }
    std::uniform_real_distribution<T> dist{0, 1. - std::numeric_limits<T>::epsilon()};
};


template<typename Particle>
auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}



template<size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static const size_t dimension = dim;

    double weight;
    double charge;

    std::array<int, dim> iCell    = ConstArray<int, dim>();
    std::array<double, dim> delta = ConstArray<double, dim>();
    std::array<double, 3> v       = ConstArray<double, 3>();


    bool operator==(Particle<dim> const& that) const
    {
        return (this->weight == that.weight) && //
               (this->charge == that.charge) && //
               (this->iCell == that.iCell) &&   //
               (this->delta == that.delta) &&   //
               (this->v == that.v);
    }

    template<std::size_t dimension>
    friend std::ostream& operator<<(std::ostream& out, const Particle<dimension>& particle);
};

template<std::size_t dim>
std::ostream& operator<<(std::ostream& out, Particle<dim> const& particle)
{
    out << "iCell(";
    for (auto c : particle.iCell)
        out << c << ",";
    out << "), delta(";
    for (auto d : particle.delta)
        out << d << ",";
    out << "), v(";
    for (auto v : particle.v)
        out << v << ",";
    out << "), charge : " << particle.charge << ", weight : " << particle.weight << '\n';
    return out;
}



} // namespace PHARE::core




#endif
