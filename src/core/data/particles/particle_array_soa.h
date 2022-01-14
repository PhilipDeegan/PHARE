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


template<std::size_t dim, std::size_t size_>
struct SoAArray
{
    std::array<double, size_> weight, charge;
    std::array<std::array<int, dim>, size_> iCell_;
    std::array<std::array<double, dim>, size_> delta;
    std::array<std::array<double, 3>, size_> v;


    auto constexpr static size() { return size_; }
};

template<std::size_t dim>
struct SoAVector
{
    using container_type = std::vector<Particle<dim>>;
    using value_type     = typename container_type::value_type;

    SoAVector() {}

    SoAVector(std::size_t size)
        : weight(size)
        , charge(size)
        , iCell_(size)
        , delta(size)
        , v(size)
    {
    }

    template<typename Particle_t>
    SoAVector(std::size_t size, Particle_t&& from)
        : weight(size, from.weight())
        , charge(size, from.charge())
        , iCell_(size, from.iCell())
        , delta(size, from.delta())
        , v(size, from.v())
    {
    }

    auto size() const { return weight.size(); }

    std::vector<double> weight, charge;
    std::vector<std::array<int, dim>> iCell_;
    std::vector<std::array<double, dim>> delta;
    std::vector<std::array<double, 3>> v;
};

template<std::size_t dim, bool _const_ = false>
struct SoASpan
{
    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<const__>>
    SoASpan(Container_int const&& _iCell, Container_double const&& _delta,
            Container_double const&& _weight, Container_double const&& _charge,
            Container_double const&& _v)
        : size_{_weight.size()}
        , weight{reinterpret_cast<double const*>(_weight.data())}
        , charge{reinterpret_cast<double const*>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim> const*>(_iCell.data())}
        , delta{reinterpret_cast<std::array<double, dim> const*>(_delta.data())}
        , v{reinterpret_cast<std::array<double, 3> const*>(_v.data())}
    {
    }
    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<!const__>>
    SoASpan(Container_int&& _iCell, Container_double&& _delta, Container_double&& _weight,
            Container_double&& _charge, Container_double&& _v)
        : size_{_weight.size()}
        , weight{reinterpret_cast<double* const>(_weight.data())}
        , charge{reinterpret_cast<double* const>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim>* const>(_iCell.data())}
        , delta{reinterpret_cast<std::array<double, dim>* const>(_delta.data())}
        , v{reinterpret_cast<std::array<double, 3>* const>(_v.data())}
    {
    }

    auto& size() const { return size_; }

    std::size_t size_;

    template<typename T>
    using container_t = std::conditional_t<_const_, T const*, T* const>;

    container_t<double> weight, charge;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta;
    container_t<std::array<double, 3>> v;
};



template<std::size_t dim, typename Super_ = SoAVector<dim>>
struct SoAParticles : public Super_
{
    using This  = SoAParticles<dim, Super_>;
    using Super = Super_;
    using Super::charge;
    using Super::delta;
    using Super::iCell_;
    using Super::size;
    using Super::v;
    using Super::weight;

    template<typename T>
    auto static constexpr is_vector()
    {
        return std::is_same_v<T, SoAVector<dim>>;
    }
    template<typename T>
    auto static constexpr is_span()
    {
        return std::is_same_v<T, SoASpan<dim, true>> or std::is_same_v<T, SoASpan<dim, false>>;
    }

    SoAParticles() {}

    template<typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    SoAParticles(std::size_t size = 0)
        : Super{size}
    {
    }

    template<typename Particle_t, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    SoAParticles(std::size_t size, Particle_t&& particle)
        : Super{size, std::forward<Particle_t>(particle)}
    {
    }

    template<typename... Args, typename S = Super, typename = std::enable_if_t<is_span<S>()>>
    SoAParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }
    template<typename... Args, typename S = Super, typename = std::enable_if_t<is_span<S>()>>
    SoAParticles(Args const&... args)
        : Super{args...}
    {
    }

    auto constexpr static size_of()
    {
        return sizeof(typename decltype(iCell_)::value_type)
               + sizeof(typename decltype(delta)::value_type)
               + sizeof(typename decltype(weight)::value_type)
               + sizeof(typename decltype(charge)::value_type)
               + sizeof(typename decltype(v)::value_type);
    }

    auto mem_size() const { return size() * size_of(); }

    auto& iCell(std::size_t i) const { return iCell_[i]; }
    auto& iCell(std::size_t i) { return iCell_[i]; }

    auto operator[](std::size_t i) { return std::make_pair(this, i); }
    auto operator[](std::size_t i) const { return std::make_pair(this, i); }

    template<typename T>
    struct iterator_t; //

    using iterator       = iterator_t<This>;
    using const_iterator = iterator_t<This const>;


    auto begin() { return iterator(*this); }
    auto begin() const { return const_iterator(*this); }
    auto end() { return iterator(*this) + size(); }
    auto end() const { return const_iterator(*this) + size(); }

    template<typename S = Super, std::enable_if_t<is_vector<S>(), int> = 0>
    void push_back(iterator const& it)
    {
        this->weight.push_back(it.weight());
        this->charge.push_back(it.charge());
        this->iCell_.push_back(it.iCell());
        this->delta.push_back(it.delta());
        this->v.push_back(it.v());
    }

    template<typename This_, typename Size, //
             typename S = Super, std::enable_if_t<is_vector<S>(), int> = 0>
    void push_back(std::pair<This_, Size> const& particle)
    {
        auto [that, idx] = particle;
        this->weight.push_back(that->weight[idx]);
        this->charge.push_back(that->charge[idx]);
        this->iCell_.push_back(that->iCell_[idx]);
        this->delta.push_back(that->delta[idx]);
        this->v.push_back(that->v[idx]);
    }
};


template<std::size_t dim, typename Super_>
template<typename T>
struct SoAParticles<dim, Super_>::iterator_t
{
    static constexpr auto dimension     = dim;
    auto static constexpr is_const      = std::is_const_v<T>;
    auto static constexpr is_contiguous = true;

    using outer_type      = std::decay_t<T>;
    using difference_type = std::size_t;
    // using value_type        = ParticleView<dim>;
    // using pointer           = ParticleView<dim>*;
    // using reference         = ParticleView<dim>&;
    using iterator_category = std::output_iterator_tag;

    iterator_t(T& particles_)
        : particles{particles_}
    {
    }

    iterator_t(iterator_t&& that)      = default;
    iterator_t(iterator_t const& that) = default;

    auto& operator++()
    {
        ++curr_pos;
        return *this;
    }
    auto& operator+=(std::size_t i)
    {
        curr_pos += i;
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

    auto& operator=(iterator_t const& that)
    {
        curr_pos = that.curr_pos;
        return *this;
    }
    auto& operator=(iterator_t&& that)
    {
        curr_pos = that.curr_pos;
        return *this;
    }

    auto operator==(iterator_t const& that) const { return curr_pos == that.curr_pos; }
    auto operator!=(iterator_t const& that) const { return curr_pos != that.curr_pos; }
    auto operator<(iterator_t const& that) const { return curr_pos < that.curr_pos; }
    // auto& operator*() { return particles.view(curr_pos); }
    // auto operator*() const { return particles.view(curr_pos); }
    auto& operator()() const { return curr_pos; }

    auto& container() { return particles; }
    auto& container() const { return particles; }


    auto& weight() { return particles.weight[curr_pos]; }
    auto& weight() const { return particles.weight[curr_pos]; }
    auto& charge() { return particles.charge[curr_pos]; }
    auto& charge() const { return particles.charge[curr_pos]; }
    auto& iCell() { return particles.iCell_[curr_pos]; }
    auto& iCell() const { return particles.iCell_[curr_pos]; }
    auto& delta() { return particles.delta[curr_pos]; }
    auto& delta() const { return particles.delta[curr_pos]; }
    auto& v() { return particles.v[curr_pos]; }
    auto& v() const { return particles.v[curr_pos]; }

    template<typename Weight>
    void weight(Weight weight)
    {
        particles.weight[curr_pos] = weight;
    }

    template<typename Charge>
    void charge(Charge charge)
    {
        particles.charge[curr_pos] = charge;
    }

    template<typename ICell>
    void iCell(ICell iCell)
    {
        particles.iCell_[curr_pos] = iCell;
    }

    template<typename Delta>
    void delta(Delta delta)
    {
        particles.delta[curr_pos] = delta;
    }

    template<typename V>
    void v(V v)
    {
        particles.v[curr_pos] = v;
    }

    T& particles; // might be const
    std::size_t curr_pos = 0;
};


template<std::size_t dim>
struct ParticleArray_SOA : public SoAParticles<dim, SoAVector<dim>>
{
    static constexpr bool is_contiguous    = true;
    static constexpr std::size_t dimension = dim;

    using This       = ParticleArray_SOA<dim>;
    using Super      = SoAParticles<dim, SoAVector<dim>>;
    using Particle_t = std::pair<Super*, std::size_t>;
    using value_type = Particle_t;
    // using view_t     = std::conditional_t<_const_, ParticleView_const<dim>, ParticleView<dim>>;

    using Super::charge;
    using Super::delta;
    using Super::iCell_;
    using Super::size;
    using Super::v;
    using Super::weight;

    // using Super::begin;
    // using Super::end;

    using iterator       = typename Super::template iterator_t<This>;
    using const_iterator = typename Super::template iterator_t<This const>;

    auto begin() { return iterator(*this); }
    auto begin() const { return const_iterator(*this); }
    auto end() { return iterator(*this) + size(); }
    auto end() const { return const_iterator(*this) + size(); }

    // using iterator       = typename Super::iterator;
    // using const_iterator = typename Super::const_iterator;

    template<std::size_t size>
    using array_type = SoAParticles<dim, SoAArray<dim, size>>;

    ParticleArray_SOA(std::size_t size = 0)
        : Super{size}
    {
    }

    // impl in particle_packer.h
    inline ParticleArray_SOA(ParticleArray<dim> const& particles);


    ParticleArray_SOA(std::size_t s, Particle<dim> from)
        : Super{s, std::forward<Particle<dim>>(from)}
    {
        // make_views();
    }

    // template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    // ParticleArray_SOA(ParticleArray_SOA const& that)
    //     : size_{that.size_}
    //     , weight{std::move(that.weight)}
    //     , charge{std::move(that.charge)}
    //     , iCell_{std::move(that.iCell_)}
    //     , delta{std::move(that.delta)}
    //     , v{std::move(that.v)}
    //     , views{}
    // {
    //     // make_views();
    // }

    // template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    // ParticleArray_SOA(ParticleArray_SOA&& that)
    //     : size_{that.size_}
    //     , weight{std::move(that.weight)}
    //     , charge{std::move(that.charge)}
    //     , iCell_{std::move(that.iCell_)}
    //     , delta{std::move(that.delta)}
    //     , v{std::move(that.v)}
    //     , views{}
    // {
    //     // make_views();
    // }




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
    // auto& view(std::size_t i) { return views[i]; }
    auto view(std::size_t i) const { return _to<ParticleView_const<dim>>(i); }

    // auto& operator[](std::size_t i) { return view(i); }
    // auto operator[](std::size_t i) const { return view(i); }

    // auto super() { return static_cast<Super*>(this); }
    // auto super() const { return const_cast<Super*>(static_cast<Super const*>(this)); }

    auto operator[](std::size_t i) { return std::make_pair(this, i); }
    auto operator[](std::size_t i) const { return std::make_pair(this, i); }

    // void push_back(iterator const& it)
    // {
    //     this->weight.push_back(it.weight());
    //     this->charge.push_back(it.charge());
    //     this->iCell_.push_back(it.iCell_());
    //     this->delta.push_back(it.delta());
    //     this->v.push_back(it.v());
    // }


    void push_back(This const& that, std::size_t idx)
    {
        this->weight.push_back(that.weight[idx]);
        this->charge.push_back(that.charge[idx]);
        this->iCell_.push_back(that.iCell_[idx]);
        this->delta.push_back(that.delta[idx]);
        this->v.push_back(that.v[idx]);
    }


    void push_back(Particle<dim> const& particle)
    { // not sure this should exist
        this->weight.push_back(particle.weight_);
        this->charge.push_back(particle.charge_);
        this->iCell_.push_back(particle.iCell_);
        this->delta.push_back(particle.delta_);
        this->v.push_back(particle.v_);
    }


    template<typename This_, typename Size>
    void push_back(std::pair<This_, Size> const& particle)
    {
        auto [that, idx] = particle;
        this->weight.push_back(that->weight[idx]);
        this->charge.push_back(that->charge[idx]);
        this->iCell_.push_back(that->iCell_[idx]);
        this->delta.push_back(that->delta[idx]);
        this->v.push_back(that->v[idx]);
    }

    template<typename... Args>
    Particle_t emplace_back(Args const&... args)
    {
        auto arg_tuple  = std::forward_as_tuple(args...);
        auto this_tuple = as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto ic) {
            auto constexpr i = ic();
            std::get<i>(this_tuple).emplace_back(std::get<i>(arg_tuple));
        });

        return (*this)[this->size() - 1];
    }


    auto as_tuple() { return std::forward_as_tuple(weight, charge, iCell_, delta, v); }
    auto as_tuple() const { return std::forward_as_tuple(weight, charge, iCell_, delta, v); }

    void clear()
    {
        std::apply([](auto&... container) { ((container.clear()), ...); }, as_tuple());
    }
    void resize(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.resize(size)), ...); }, as_tuple());
    }
    void reserve(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.reserve(size)), ...); }, as_tuple());
    }
    void swap(This& that)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = that.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto ic) {
            auto constexpr i = ic();
            std::get<i>(this_tuple).swap(std::get<i>(that_tuple));
        });
    }


    template<typename V>
    void swap_(V& a, V& b)
    {
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


    ParticleArray_SOA(iterator start, iterator end);

    iterator erase(iterator first, iterator last)
    {
        throw std::runtime_error("fix me");
        // clean_ = false;
        // // return particles_.erase(first, last);
        // return iterator{particles_.erase(first, last), this};
    }

    struct Sorter;
    void sort() { Sorter{*this}(); }
    void sort(std::int64_t s, std::int64_t e) { Sorter{*this}(s, e); }
};



template<std::size_t dim>
struct ParticleArray_SOA<dim>::Sorter
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

    ParticleArray_SOA<dim>& self;
};


template<std::size_t dim>
using ParticleArray_SOA_t = ParticleArray_SOA<dim>;

template<std::size_t dim, bool _const_ = false>
using ParticleArray_SOAView = SoAParticles<dim, SoASpan<dim, _const_>>;



template<std::size_t dim>
void empty(ParticleArray_SOA<dim>& array)
{
    array.clear();
}

template<std::size_t dim>
void swap(ParticleArray_SOA<dim>& array1, ParticleArray_SOA<dim>& array2)
{
    array1.swap(array2);
}

} // namespace PHARE::core


namespace std
{
using namespace PHARE::core;

template<std::size_t dim, typename Super>
Particle<dim> copy(std::pair<SoAParticles<dim, Super>*, std::size_t> const& pair)
{
    auto const& [from, idx] = pair;
    return {from.weight[idx], from.charge[idx], from.iCell[idx], from.delta[idx], from.v[idx]};
}

template<std::size_t dim>
Particle<dim> copy(std::pair<ParticleArray_SOA<dim>*, std::size_t> const& pair)
{
    auto const& [from, idx] = pair;
    return {from.weight[idx], from.charge[idx], from.iCell[idx], from.delta[idx], from.v[idx]};
}



template<template<typename> typename Iterator, typename O, typename Inserter>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void> copy(Iterator<O> const& a,
                                                                 Iterator<O> const& b, Inserter i)
{
    throw std::runtime_error("Fix me");
}

template<template<typename> typename Iterator, typename O, typename Inserter, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void>
copy_if(Iterator<O> src_begin, Iterator<O> src_end, Inserter i, Fn check)
{
    auto& particles = src_begin.particles;

    for (; src_begin != src_end; ++src_begin)
        if (check(src_begin))
        {
            // auto const& [_, idx] = particles[src_begin.curr_pos];
            // auto pair            = std::make_pair(&particles, idx);
            *i = particles[src_begin.curr_pos];
            ++i;
        }
}


template<template<typename> typename Iterator, typename O, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, Iterator<O>>
partition(Iterator<O> a, Iterator<O> b, Fn check)
{
    throw std::runtime_error("Fix me");
    return a;
}

template<template<typename> typename Iterator, typename O, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void>
transform(Iterator<O> src_begin, Iterator<O> src_end, Iterator<O> dst_begin, Fn transform)
{
    if (&src_begin.particles != &dst_begin.particles)
        throw std::runtime_error("no");

    for (; src_begin != src_end; ++src_begin, ++dst_begin)
        /*dst_begin = */ transform(src_begin); // transform edits in place
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, std::size_t> distance(Iterator<O> const& a,
                                                                            Iterator<O> const& b)
{
    assert(&a.particles == &b.particles);
    assert(a.curr_pos <= b.curr_pos);

    return b.curr_pos - a.curr_pos;
}

template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, std::size_t> operator-(Iterator<O> const& a,
                                                                             Iterator<O> const& b)
{
    assert(&a.particles == &b.particles);
    assert(b.curr_pos <= a.curr_pos);

    return a.curr_pos - b.curr_pos;
}


template<std::size_t dim>
void sort(ParticleArray_SOA<dim>& particles)
{
    particles.sort();
}


} // namespace std


namespace PHARE::core
{
template<std::size_t dim>
// template<bool OS /*enable_if_t*/>
ParticleArray_SOA<dim>::ParticleArray_SOA(ParticleArray_SOA<dim>::iterator start,
                                          ParticleArray_SOA<dim>::iterator end)
    : Super{std::distance(start, end)}
// , weight(size_)
// , charge(size_)
// , iCell_(size_)
// , delta(size_)
// , v(size_)
// , views{}
{
    // // make_views();
    throw std::runtime_error("fix me");
}

} // namespace PHARE::core


#endif