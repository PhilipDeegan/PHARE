#ifndef PHARE_CORE_DATA_PARTICLES_CONTIGUOUS_H
#define PHARE_CORE_DATA_PARTICLES_CONTIGUOUS_H

#include "particle.h"
#include "particle_array.h"

namespace PHARE::core
{
template<typename Float, std::size_t dim>
struct ContiguousParticle
{
    int* iCell    = nullptr;
    float* delta  = nullptr;
    Float *weight = nullptr, *charge = nullptr, *v = nullptr;
};



template<typename Float, std::size_t dim, bool OwnedState = true>
struct ContiguousParticles
{
    static constexpr bool is_contiguous    = true;
    static constexpr std::size_t dimension = dim;
    using ContiguousParticles_             = ContiguousParticles<Float, dim, OwnedState>;

    template<typename T>
    using container_t = std::conditional_t<OwnedState, std::vector<T>, Span<T>>;

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ContiguousParticles(std::size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
    {
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    inline ContiguousParticles(ParticleArray<Float, dim> const& particles);

    template<typename Container_int, typename Container_float, typename Container_double>
    ContiguousParticles(Container_int&& _iCell, Container_float&& _delta,
                        Container_double&& _weight, Container_double&& _charge,
                        Container_double&& _v)
        : iCell{_iCell}
        , delta{_delta}
        , weight{_weight}
        , charge{_charge}
        , v{_v}
    {
    }

    std::size_t size() const { return weight.size(); }

    template<std::size_t S, typename T>
    static std::array<T, S>* _array_cast(T const* array)
    {
        return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
    }

    template<typename Return>
    Return _to(std::size_t i)
    {
        return {
            *const_cast<Float*>(weight.data() + i),      //
            *const_cast<Float*>(charge.data() + i),      //
            *_array_cast<dim>(iCell.data() + (dim * i)), //
            *_array_cast<dim>(delta.data() + (dim * i)), //
            *_array_cast<3>(v.data() + (3 * i)),
        };
    }

    auto copy(std::size_t i) { return _to<Particle<Float, dim>>(i); }
    auto view(std::size_t i) { return _to<ParticleView<Float, dim>>(i); }

    auto operator[](std::size_t i) const { return view(i); }
    auto operator[](std::size_t i) { return view(i); }

    struct iterator
    {
        iterator(ContiguousParticles_* particles)
        {
            for (std::size_t i = 0; i < particles->size(); i++)
                views.emplace_back((*particles)[i]);
        }

        iterator& operator++()
        {
            ++curr_pos;
            return *this;
        }

        bool operator!=(iterator const& other) const { return curr_pos != views.size(); }
        auto& operator*() { return views[curr_pos]; }
        auto& operator*() const { return views[curr_pos]; }

        std::size_t curr_pos = 0;
        std::vector<ParticleView<Float, dim>> views;
    };

    auto begin() { return iterator(this); }
    auto cbegin() const { return iterator(this); }

    auto end() { return iterator(this); }
    auto cend() const { return iterator(this); }

    container_t<int> iCell;
    container_t<float> delta;
    container_t<Float> weight, charge, v;
};

template<typename Float, std::size_t dim>
using ContiguousParticlesView = ContiguousParticles<Float, dim, /*OwnedState=*/false>;


template<typename Float, std::size_t dim>
struct EMParticleView : ContiguousParticle<Float, dim>
{
    Float *E = nullptr, *B = nullptr;
};

template<typename Float, std::size_t dim>
struct EMContiguousParticles : ContiguousParticles<Float, dim>
{
    using Super                        = ContiguousParticles<Float, dim>;
    static constexpr std::size_t THREE = 3;

    EMContiguousParticles(std::size_t s)
        : Super(s)
        , E(s * THREE)
        , B(s * THREE)
    {
    }

    inline EMContiguousParticles(ParticleArray<Float, dim> const& particles);

    void swap(std::size_t idx0, std::size_t idx1)
    {
        Super::swap(idx0, idx1);
        Super::_swap(idx0, idx1, E, THREE);
        Super::_swap(idx0, idx1, B, THREE);
    }
    EMParticleView<Float, dim> particle(std::size_t i) const
    {
        return {
            Super::iCell.data() + (dim * i), //
            Super::delta.data() + (dim * i), //
            Super::weight.data() + i,        //
            Super::charge.data() + i,        //
            Super::v.data() + (THREE * i),   //
            E.data() + (THREE * i),          //
            B.data() + (THREE * i),
        };
    }

    std::vector<Float> E, B;
};

template<typename Float, std::size_t dim>
class ParticlePacker
{
    static constexpr std::size_t THREE = 3;

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


    void _pack(ContiguousParticles<Float, dim>& to, std::size_t idx)
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
        for (std::size_t idx = 0; idx < particles_.size(); idx++)
            _pack(to, idx);
    }

    void pack(EMContiguousParticles<Float, dim>& to)
    {
        for (std::size_t idx = 0; idx < particles_.size(); idx++)
        {
            _pack(to, idx);
            _copy(&particles_[idx].Ex, idx, THREE, to.E);
            _copy(&particles_[idx].Bx, idx, THREE, to.B);
        }
    }

private:
    ParticleArray<Float, dim> const& particles_;

    template<typename Array, typename Vector>
    void _copy(Array* from, std::size_t idx, std::size_t size, Vector& to)
    {
        std::copy(from, from + size, to.begin() + (idx * size));
    }
};


template<typename Float, std::size_t dim, bool OwnedState>
template<bool OS, typename>
ContiguousParticles<Float, dim, OwnedState>::ContiguousParticles(
    ParticleArray<Float, dim> const& particles)
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
