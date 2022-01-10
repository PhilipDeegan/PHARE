#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_H

#include "core/data/particles/particle.h"

namespace PHARE::core
{
// forward decl
template<std::size_t dim>
class ParticleArray;



template<std::size_t dim, bool _const_ = false>
struct ParticleViewBase
{
    static constexpr std::size_t dimension = dim;
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");

    using This = ParticleViewBase<dim, _const_>;

    template<typename T>
    using if_const_t = std::conditional_t<_const_, T const, T>;

    ParticleViewBase(if_const_t<double>& weight_, if_const_t<double>& charge_,
                     if_const_t<std::array<int, dim>>& iCell_,
                     if_const_t<std::array<double, dim>>& delta_,
                     if_const_t<std::array<double, 3>>& v_)
        : weight{weight_}
        , charge{charge_}
        , iCell{iCell_}
        , delta{delta_}
        , v{v_}
    {
    }

    ParticleViewBase(ParticleViewBase const& that)
        : weight{that.weight}
        , charge{that.charge}
        , iCell{that.iCell}
        , delta{that.delta}
        , v{that.v}
    {
    }

    auto operator<(This const& that) const { return iCell < that.iCell; }

    if_const_t<double>& weight;
    if_const_t<double>& charge;
    if_const_t<std::array<int, dim>>& iCell;
    if_const_t<std::array<double, dim>>& delta;
    if_const_t<std::array<double, 3>>& v;
};

template<std::size_t dim>
struct ParticleView : ParticleViewBase<dim>
{
};

template<std::size_t dim>
struct ParticleView_const : ParticleViewBase<dim, true>
{
};

template<std::size_t dim, typename T>
inline constexpr auto is_phare_particle_view_type =   //
    std::is_same_v<ParticleViewBase<dim>, T> or       //
    std::is_same_v<ParticleViewBase<dim, true>, T> or //
    std::is_same_v<ParticleView<dim>, T> or           //
    std::is_same_v<ParticleView_const<dim>, T>;


template<std::size_t dim, bool OwnedState = true, bool _const_ = false>
struct ParticleArray_SOA
{
    using This       = ParticleArray_SOA<dim, OwnedState>;
    using value_type = ParticleView<dim>;
    using view_t     = std::conditional_t<_const_, ParticleView_const<dim>, ParticleView<dim>>;

    static constexpr bool is_contiguous    = true;
    static constexpr std::size_t dimension = dim;

    auto make_views()
    {
        views = generate([&](auto i) { return this->_to<view_t>(i); }, size_);
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA(std::size_t s)
        : size_{s}
        , weight(s)
        , charge(s)
        , iCell_(s)
        , delta(s)
        , v(s)
    {
        make_views();
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    // impl in particle_packer.h
    inline ParticleArray_SOA(ParticleArray<dim> const& particles);


    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA(std::size_t s, Particle<dim> from)
        : size_{s}
        , weight(s, from.weight)
        , charge(s, from.charge)
        , iCell_(s, from.iCell)
        , delta(s, from.delta)
        , v(s, from.v)
    {
        make_views();
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA(ParticleArray_SOA const& that)
        : size_{that.size_}
        , weight{std::move(that.weight)}
        , charge{std::move(that.charge)}
        , iCell_{std::move(that.iCell_)}
        , delta{std::move(that.delta)}
        , v{std::move(that.v)}
        , views{}
    {
        make_views();
    }

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ParticleArray_SOA(ParticleArray_SOA&& that)
        : size_{that.size_}
        , weight{std::move(that.weight)}
        , charge{std::move(that.charge)}
        , iCell_{std::move(that.iCell_)}
        , delta{std::move(that.delta)}
        , v{std::move(that.v)}
        , views{}
    {
        make_views();
    }


    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<const__>>
    ParticleArray_SOA(Container_int const&& _iCell, Container_double const&& _delta,
                      Container_double const&& _weight, Container_double const&& _charge,
                      Container_double const&& _v)
        : size_{_weight.size()}
        , weight{reinterpret_cast<double const*>(_weight.data())}
        , charge{reinterpret_cast<double const*>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim> const*>(_iCell.data())}
        , delta{reinterpret_cast<std::array<double, dim> const*>(_delta.data())}
        , v{reinterpret_cast<std::array<double, 3> const*>(_v.data())}
    {
        make_views();
    }
    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<!const__>>
    ParticleArray_SOA(Container_int&& _iCell, Container_double&& _delta, Container_double&& _weight,
                      Container_double&& _charge, Container_double&& _v)
        : size_{_weight.size()}
        , weight{reinterpret_cast<double* const>(_weight.data())}
        , charge{reinterpret_cast<double* const>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim>* const>(_iCell.data())}
        , delta{reinterpret_cast<std::array<double, dim>* const>(_delta.data())}
        , v{reinterpret_cast<std::array<double, 3>* const>(_v.data())}
    {
        make_views();
    }

    auto mem_size() const
    {
        return size_
               * (sizeof(typename decltype(iCell_)::value_type) //
                  + sizeof(typename decltype(delta)::value_type)
                  + sizeof(typename decltype(weight)::value_type)
                  + sizeof(typename decltype(charge)::value_type)
                  + sizeof(typename decltype(v)::value_type));
    }

    template<typename S>
    bool _check(S i) const // templated so doesn't exist in release mode
    {
        return i < size();
    }


    auto& iCell(std::size_t i) const
    {
        assert(_check(i));
        return iCell_[i];
    }
    auto& iCell(std::size_t i)
    {
        assert(_check(i));
        return iCell_[i];
    }

    template<typename Return>
    Return _to(std::size_t i)
    {
        return {{weight[i], charge[i], iCell_[i], delta[i], v[i]}};
    }
    template<typename Return>
    Return _to(std::size_t i) const
    {
        return {{weight[i], charge[i], iCell_[i], delta[i], v[i]}};
    }

    auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
    auto& view(std::size_t i) { return views[i]; }
    auto view(std::size_t i) const { return _to<ParticleView_const<dim>>(i); }

    auto& operator[](std::size_t i) { return view(i); }
    auto operator[](std::size_t i) const { return view(i); }


    template<typename T>
    // this is to allow interop with SoA, but is to be avoided as it's not effecient
    struct iterator_impl
    {
    public:
        auto constexpr static is_const = std::is_const_v<T>;

        using outer_type        = This;
        using difference_type   = std::size_t;
        using value_type        = ParticleView<dim>;
        using pointer           = ParticleView<dim>*;
        using reference         = ParticleView<dim>&;
        using iterator_category = std::output_iterator_tag;

        iterator_impl(T& particles_)
            : particles{particles_}
        {
        }

        iterator_impl(iterator_impl&& that)      = default;
        iterator_impl(iterator_impl const& that) = default;

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

        auto& operator=(iterator_impl const& that)
        {
            curr_pos = that.curr_pos;
            return *this;
        }
        auto& operator=(iterator_impl&& that)
        {
            curr_pos = that.curr_pos;
            return *this;
        }

        auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
        auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
        auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
        auto& operator*() { return particles[curr_pos]; }
        auto operator*() const { return particles[curr_pos]; }
        auto& operator()() const { return curr_pos; }

        auto& container() { return particles; }
        auto& container() const { return particles; }

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


    template<typename V>
    void swap_(V& a, V& b)
    {
        // static thread_local V t;
        // t = a;
        // a = b;
        // b = t;

        std::swap(a, b);
    }

    void swap(std::size_t const& a, std::size_t const& b)
    {
        // if (a == b)
        //     return;
        swap_(weight[a], weight[b]);
        swap_(charge[a], charge[b]);
        swap_(iCell_[a], iCell_[b]);
        swap_(delta[a], delta[b]);
        swap_(v[a], v[b]);
    }

    struct Sorter
    {
        template<typename V>
        auto el_wise_less(V const& v0, V const& v1)
        {
            for (std::int16_t i = 0; i < v0.size(); ++i)
                if (v0[i] < v1[i])
                    return true;
                else if (v0[i] != v1[i])
                    return false;
            return false;
        }
        template<typename V>
        auto el_wise_gr8r(V const& v0, V const& v1)
        {
            for (std::int16_t i = 0; i < v0.size(); ++i)
                if (v0[i] > v1[i])
                    return true;
                else if (v0[i] != v1[i])
                    return false;
            return false;
        }

        void operator()(std::int64_t l, std::int64_t r)
        {
            auto i = l;
            auto j = r;

            auto const& iCells = self.iCell_;
            auto const half    = iCells[(l + r) / 2];

            do
            {
                while (el_wise_less(iCells[i], half))
                    i++;
                while (el_wise_gr8r(iCells[j], half))
                    j--;

                if (i <= j)
                {
                    self.swap(i, j);

                    i++;
                    j--;
                }
            } while (i <= j);

            if (l < j)
                (*this)(l, j);
            if (i < r)
                (*this)(i, r);
        }

        void operator()() { (*this)(0, self.size() - 1); }

        This& self;
    };

    void sort() { Sorter{*this}(); }
    void sort(std::int64_t s, std::int64_t e) { Sorter{*this}(s, e); }


    template<typename T>
    using container_t = std::conditional_t<OwnedState, std::vector<T>,
                                           std::conditional_t<_const_, T const*, T* const>>;


    std::size_t size_;
    container_t<double> weight, charge;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta;
    container_t<std::array<double, 3>> v;

    std::vector<view_t> views;
};


template<std::size_t dim>
using ParticleArray_SOA_t = ParticleArray_SOA<dim>;

template<std::size_t dim, bool _const_ = false>
using ParticleArray_SOAView = ParticleArray_SOA<dim, /*OwnedState=*/false, _const_>;


template<std::size_t dim>
void swap(ParticleView<dim>& a, ParticleView<dim>& b)
{
    std::swap(a.weight, b.weight);
    std::swap(a.charge, b.charge);
    std::swap(a.iCell, b.iCell);
    std::swap(a.delta, b.delta);
    std::swap(a.v, b.v);
}


} // namespace PHARE::core


namespace std
{
using namespace PHARE::core;

template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::outer_type::is_contiguous, std::size_t>
distance(Iterator<O> const& a, Iterator<O> const& b)
{
    assert(&a.particles == &b.particles);
    assert(a.curr_pos < b.curr_pos);

    return b.curr_pos - a.curr_pos;
}

template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::outer_type::is_contiguous, std::size_t>
operator-(Iterator<O> const& a, Iterator<O> const& b)
{
    assert(&a.particles == &b.particles);
    assert(b.curr_pos <= a.curr_pos);

    return a.curr_pos - b.curr_pos;
}


template<std::size_t dim, bool OwnedState = true>
void sort(ParticleArray_SOA<dim, OwnedState>& particles)
{
    particles.sort();
}


} // namespace std



#endif