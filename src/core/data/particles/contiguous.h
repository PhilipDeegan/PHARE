#ifndef PHARE_CORE_DATA_PARTICLES_CONTIGUOUS_H
#define PHARE_CORE_DATA_PARTICLES_CONTIGUOUS_H

#include "particle_array.h"

namespace PHARE::core
{
template<typename Float, size_t dim>
struct ContiguousParticle
{
    int* iCell    = nullptr;
    float* delta  = nullptr;
    Float *weight = nullptr, *charge = nullptr, *v = nullptr;
};


template<typename Float, std::size_t dim>
struct ContiguousParticles
{
    static constexpr size_t ONE   = 1;
    static constexpr size_t THREE = 3;

    ContiguousParticles(size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * THREE)
        , size_(s)
    {
    }
    inline ContiguousParticles(ParticleArray<Float, dim> const& particles);

    auto size() const { return size_; }

    void swap(size_t idx0, size_t idx1)
    {
        _swap(idx0, idx1, iCell, dim);
        _swap(idx0, idx1, delta, dim);
        _swap(idx0, idx1, v, THREE);
        _swap(idx0, idx1, weight, ONE);
        _swap(idx0, idx1, charge, ONE);
    }

    template<typename T>
    void _swap(size_t idx0, size_t idx1, T*& v, size_t size)
    {
        std::swap_ranges(v + idx0, v + idx0 + size, v + idx1);
    }

    template<typename T>
    void _swap(size_t idx0, size_t idx1, std::vector<T>& v, size_t size)
    {
        _swap(idx0, idx1, v.data(), size);
    }

    auto particle(size_t i) const
    {
        return ContiguousParticle<Float, dim>{
            iCell.data() + (dim * i), //
            delta.data() + (dim * i), //
            weight.data() + i,        //
            charge.data() + i,        //
            v.data() + (THREE * i),
        };
    }

    std::vector<int> iCell;
    std::vector<float> delta;
    std::vector<Float> weight, charge, v;
    size_t size_;
};

template<typename Float, size_t dim>
struct EMContiguousParticle : ContiguousParticle<Float, dim>
{
    double *E = nullptr, *B = nullptr;
};

template<typename Float, std::size_t dim>
struct EMContiguousParticles : ContiguousParticles<Float, dim>
{
    using Super                   = ContiguousParticles<Float, dim>;
    static constexpr size_t THREE = Super::THREE;

    EMContiguousParticles(size_t s)
        : Super(s)
        , E(s * THREE)
        , B(s * THREE)
    {
    }

    inline EMContiguousParticles(ParticleArray<Float, dim> const& particles);

    void swap(size_t idx0, size_t idx1)
    {
        Super::swap(idx0, idx1);
        Super::_swap(idx0, idx1, E, THREE);
        Super::_swap(idx0, idx1, B, THREE);
    }
    auto particle(size_t i) const
    {
        return EMContiguousParticle<Float, dim>{
            Super::iCell.data() + (dim * i), //
            Super::delta.data() + (dim * i), //
            Super::weight.data() + i,        //
            Super::charge.data() + i,        //
            Super::v.data() + (THREE * i),   //
            E.data() + (THREE * i),          //
            B.data() + (THREE * i),
        };
    }

    std::vector<double> E, B;
};

template<typename Float, size_t dim>
class ParticlePacker
{
    static constexpr size_t THREE = ContiguousParticles<Float, dim>::THREE;

public:
    ParticlePacker(ParticleArray<Float, dim> const& particles)
        : particles_{particles}
    {
    }

    static auto keys()
    {
        return std::vector<std::string>{"weight", "charge", "iCell", "delta", "v"};
    }

    static auto forward_as_tuple(Particle<Float, dim> const& part)
    {
        return std::forward_as_tuple(part.weight, part.charge, part.iCell, part.delta, part.v);
    }

    static auto empty()
    {
        Particle<Float, dim> part;
        return forward_as_tuple(part);
    }


    void _pack(ContiguousParticles<Float, dim>& to, size_t idx)
    {
        auto tuple     = forward_as_tuple(particles_[idx]);
        to.weight[idx] = std::get<0>(tuple);
        to.charge[idx] = std::get<1>(tuple);
        _copy(std::get<2>(tuple).data(), idx, dim, to.iCell);
        _copy(std::get<3>(tuple).data(), idx, dim, to.delta);
        _copy(std::get<4>(tuple).data(), idx, THREE, to.v);
    }

    void pack(ContiguousParticles<Float, dim>& to)
    {
        for (size_t idx = 0; idx < particles_.size(); idx++)
            _pack(to, idx);
    }

    void pack(EMContiguousParticles<Float, dim>& to)
    {
        for (size_t idx = 0; idx < particles_.size(); idx++)
        {
            _pack(to, idx);
            _copy(&particles_[idx].Ex, idx, THREE, to.E);
            _copy(&particles_[idx].Bx, idx, THREE, to.B);
        }
    }

private:
    ParticleArray<Float, dim> const& particles_;

    template<typename Array, typename Vector>
    void _copy(Array* from, size_t idx, size_t size, Vector& to)
    {
        std::copy(from, from + size, to.begin() + (idx * size));
    }
};



template<typename Float, std::size_t dim>
ContiguousParticles<Float, dim>::ContiguousParticles(ParticleArray<Float, dim> const& particles)
    : ContiguousParticles{particles.size()}
{
    ParticlePacker<Float, dim>{particles}.pack(*this);
}

template<typename Float, std::size_t dim>
EMContiguousParticles<Float, dim>::EMContiguousParticles(ParticleArray<Float, dim> const& particles)
    : EMContiguousParticles{particles.size()}
{
    ParticlePacker<Float, dim>{particles}.pack(*this);
}


} // namespace PHARE::core

#endif
