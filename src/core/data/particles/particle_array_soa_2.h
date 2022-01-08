#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_2_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_2_H

#include "core/data/particles/particle_array_soa.h"

namespace PHARE::core
{
template<std::size_t dim, bool OwnedState = true>
struct ParticleArray_SOA_2
{
    using This       = ParticleArray_SOA_2<dim, OwnedState>;
    using value_type = ParticleView<dim>;

    static constexpr bool is_contiguous    = true;
    static constexpr std::size_t dimension = dim;

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA_2(std::size_t s)
        : size_{s}
        , iCell(s)
        , delta(s)
        , weight(s)
        , charge(s)
        , v(s)
    {
    }


    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    // impl in particle_packer.h
    inline ParticleArray_SOA_2(ParticleArray<dim> const& particles);


    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA_2(std::size_t s, Particle<dim> from)
        : size_{s}
        , weight(s, from.weight)
        , charge(s, from.charge)
        , iCell(s, from.iCell)
        , delta(s, from.delta)
        , v(s, from.v)
    {
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA_2(ParticleArray_SOA_2 const&) = delete;

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA_2(ParticleArray_SOA_2&& that) = delete;


    template<typename Container_int, typename Container_double>
    ParticleArray_SOA_2(Container_int&& _iCell, Container_double&& _delta,
                        Container_double&& _weight, Container_double&& _charge,
                        Container_double&& _v)
        : size_{_weight.size()}
        , iCell{_iCell}
        , delta{_delta}
        , weight{_weight}
        , charge{_charge}
        , v{_v}
    {
    }

    template<typename Return>
    Return _to(std::size_t i)
    {
        return {weight[i], charge[i], iCell[i], delta[i], v[i]};
    }
    template<typename Return>
    Return _to(std::size_t i) const
    {
        return {weight[i], charge[i], iCell[i], delta[i], v[i]};
    }

    auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
    auto view(std::size_t i) { return _to<Particle<dim>>(i); }
    auto view(std::size_t i) const { return _to<ParticleView<dim, true> const>(i); }

    auto operator[](std::size_t i) { return view(i); }
    auto operator[](std::size_t i) const { return view(i); }


    template<typename T>
    // this is to allow interop with SoA, but is to be avoided as it's not effecient
    struct iterator_impl
    {
    public:
        auto constexpr static is_const = std::is_const_v<T>;

        using outer_type        = This;
        using difference_type   = std::size_t;
        using value_type        = iterator_impl<T>;
        using pointer           = iterator_impl<T>*;
        using reference         = iterator_impl<T>&;
        using iterator_category = std::output_iterator_tag;

        iterator_impl(T& particles_)
            : particles{particles_}
        {
        }

        auto& operator++()
        {
            ++curr_pos;
            return *this;
        }
        auto& operator--()
        {
            --curr_pos;
            return *this;
        }
        auto operator+(std::size_t i)
        {
            auto copy = *this;
            copy.curr_pos += i;
            return copy;
        }
        auto operator-(std::size_t i)
        {
            auto copy = *this;
            copy.curr_pos -= i;
            return copy;
        }

        auto operator=(iterator_impl const& that) { return curr_pos = that.curr_pos; }
        auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
        auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
        auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
        auto& operator*() { return *this; }
        auto& operator*() const { return *this; }
        auto& operator()() const { return curr_pos; }

        auto& iCell() { return particles.iCell[curr_pos]; }
        auto& iCell() const { return particles.iCell[curr_pos]; }

        T& particles; // might be const
        std::size_t curr_pos = 0;
    };

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    auto begin() { return iterator(*this); }
    auto begin() const { return const_iterator(*this); }
    auto end() { return iterator(*this) + size(); }
    auto end() const { return const_iterator(*this) + size(); }

    auto size() const { return size_; }

    void swap(std::size_t const& a, std::size_t const& b)
    {
        std::swap(weight[a], weight[b]);
        std::swap(charge[a], charge[b]);
        std::swap(iCell[a], iCell[b]);
        std::swap(delta[a], delta[b]);
        std::swap(v[a], v[b]);
    }

    template<typename T>
    using container_t = std::conditional_t<OwnedState, std::vector<T>, Span<T>>;

    container_t<double> weight, charge;
    container_t<std::array<int, dim>> iCell;
    container_t<std::array<double, dim>> delta;
    container_t<std::array<double, 3>> v;
    std::size_t size_;

    // std::vector<ParticleView<dim, false>> views;
};


template<std::size_t dim>
using ParticleArray_SOA_2_t = ParticleArray_SOA_2<dim>;

template<std::size_t dim>
using ParticleArray_SOA_2View = ParticleArray_SOA_2<dim, /*OwnedState=*/false>;


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::outer_type::is_contiguous> swap(Iterator<O>& a,
                                                                       Iterator<O>& b)
{
    assert(&a.particles == &b.particles);

    a.particles.swap(a.curr_pos, b.curr_pos);
}
} // namespace PHARE::core


namespace std
{
using namespace PHARE::core;

// template<template<typename> typename Iterator, typename O>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::outer_type::is_contiguous, std::size_t>
// distance(Iterator<O> const& a, Iterator<O> const& b)
// {
//     assert(&a.particles == &b.particles);
//     assert(a.curr_pos < b.curr_pos);

//     return b.curr_pos - a.curr_pos;
// }

// template<template<typename> typename Iterator, typename O>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::outer_type::is_contiguous, std::size_t>
// operator-(Iterator<O> const& a, Iterator<O> const& b)
// {
//     assert(&a.particles == &b.particles);
//     assert(b.curr_pos <= a.curr_pos);

//     return a.curr_pos - b.curr_pos;
// }


template<std::size_t dim, bool OwnedState = true>
void sort(ParticleArray_SOA_2<dim, OwnedState>& particles)
{
    std::sort(particles.begin(), particles.end(), [](auto const& x0, auto const& y0) {
        auto& x = x0.iCell();
        auto& y = y0.iCell();

        // for (std::size_t i = 0; i < dim; ++i)
        //     if (x[i] < y[i])
        //         return true;
        // return false;
        return x < y;
    });
}


} // namespace std



#endif