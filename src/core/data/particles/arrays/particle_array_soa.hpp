#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP

#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"


namespace PHARE::core
{
template<std::size_t dim, std::size_t size_>
struct SoAArray
{
    auto static constexpr storage_mode = StorageMode::ARRAY;
    auto static constexpr dimension    = dim;
    auto static constexpr is_vector    = false;

    std::array<double, size_> weight_, charge_;
    std::array<std::array<int, dim>, size_> iCell_;
    std::array<std::array<double, dim>, size_> delta_;
    std::array<std::array<double, 3>, size_> v_;

    auto constexpr static size() { return size_; }
};


// used in SOA_ParticelArray::operator[](std::size_t) for std algorithm interactions
//  allows copying between vector and array impls without iterators
//   could be a better way maybe but it's not so simple to change
template<std::size_t dim>
using SoAParticle_crt = std::tuple<double const&,                  //  weight
                                   double const&,                  // charge
                                   std::array<int, dim> const&,    // iCell
                                   std::array<double, dim> const&, // delta
                                   std::array<double, 3> const&    // v

                                   >;
// template<std::size_t dim>
// using SoAParticle_rt = std::tuple<double&,                  //  weight
//                                   double&,                  // charge
//                                   std::array<int, dim>&,    // iCell
//                                   std::array<double, dim>&, // delta
//                                   std::array<double, 3>&   // v
//                                   >;



// used when the memory is owned elsewhere, e.g. numpy arrays
template<std::size_t dim, bool _const_ = false>
struct SoASpan
{
    auto static constexpr storage_mode = StorageMode::SPAN;
    auto static constexpr dimension    = dim;

    template<typename T>
    using container_t = std::conditional_t<_const_, T const*, T* const>;


    template<typename C0, typename C1>
    SoASpan(C0 _iCell, C1 _delta, C1 _weight, C1 _charge, C1 _v)
        : size_{_weight.size()}
        , weight_{reinterpret_cast<double*>(_weight.data())}
        , charge_{reinterpret_cast<double*>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim>*>(_iCell.data())}
        , delta_{reinterpret_cast<std::array<double, dim>*>(_delta.data())}
        , v_{reinterpret_cast<std::array<double, 3>*>(_v.data())}
    {
    }

    template<typename ParticleArray>
    SoASpan(ParticleArray&& array)
        : size_{array.size()}
        , weight_{reinterpret_cast<double* const>(&array.weight()[0])}
        , charge_{reinterpret_cast<double* const>(&array.charge()[0])}
        , iCell_{reinterpret_cast<std::array<int, dim>* const>(&array.iCell()[0])}
        , delta_{reinterpret_cast<std::array<double, dim>* const>(&array.delta()[0])}
        , v_{reinterpret_cast<std::array<double, 3>* const>(&array.v()[0])}
    {
    }

    auto& size() const { return size_; }
    void resize(std::size_t s)
    {
        PHARE_ASSERT(s <= size_); // can't be bigger
        size_ = s;
    }

    std::size_t size_;
    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;
};



template<std::size_t dim, auto alloc_mode_ = AllocatorMode::CPU>
struct SoAVector
{
    friend class SoASpan<dim>;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr dimension    = dim;

    // using container_type             = std::vector<Particle<dim>>;
    // using value_type                 = typename container_type::value_type;

    template<typename Type>
    using container_t = typename Vector<Type, alloc_mode>::vector_t;

    SoAVector() {}

    SoAVector(std::size_t size)
        : weight_(size)
        , charge_(size)
        , iCell_(size)
        , delta_(size)
        , v_(size)
    {
    }

    template<typename Particle_t>
    SoAVector(std::size_t size, Particle_t&& from)
        : weight_(size, from.weight())
        , charge_(size, from.charge())
        , iCell_(size, from.iCell())
        , delta_(size, from.delta())
        , v_(size, from.v())
    {
    }

    auto size() const { return weight_.size(); }


    void resize(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.resize(size)), ...); }, as_tuple());
    }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }

    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;
};


template<typename Super_>
class SoAParticles : public Super_
{
    template<typename T>
    struct iterator_impl;

public:
    using Super                        = Super_;
    using This                         = SoAParticles<Super>;
    auto static constexpr dimension    = Super_::dimension;
    auto static constexpr layout_mode  = LayoutMode::SoA;
    auto static constexpr storage_mode = Super::storage_mode;
    using Particle_t                   = SoAParticle_crt<dimension>;
    using Super::size;

    template<std::size_t size>
    using array_type = SoAParticles<SoAArray<dimension, size>>;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    // public for pybind but avoid otherwise
    using Super::charge_;
    using Super::delta_;
    using Super::iCell_;
    // using Super::resize;
    using Super::v_;
    using Super::weight_;


    template<typename... Args>
    SoAParticles(Args&&... args) _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    SoAParticles(iterator start, iterator end); // impl @ bottom of this file

    SoAParticles(This const& that)    = default;
    SoAParticles(This&& that)         = default;
    This& operator=(This&& that)      = default;
    This& operator=(This const& that) = default;




    auto& weight(std::size_t i) const _PHARE_ALL_FN_ { return weight_[i]; }
    auto& weight(std::size_t i) _PHARE_ALL_FN_ { return weight_[i]; }

    auto& charge(std::size_t i) const _PHARE_ALL_FN_ { return charge_[i]; }
    auto& charge(std::size_t i) _PHARE_ALL_FN_ { return charge_[i]; }

    auto& iCell(std::size_t i) const _PHARE_ALL_FN_ { return iCell_[i]; }
    auto& iCell(std::size_t i) _PHARE_ALL_FN_ { return iCell_[i]; }

    auto& delta(std::size_t i) const _PHARE_ALL_FN_ { return delta_[i]; }
    auto& delta(std::size_t i) _PHARE_ALL_FN_ { return delta_[i]; }

    auto& v(std::size_t i) const _PHARE_ALL_FN_ { return v_[i]; }
    auto& v(std::size_t i) _PHARE_ALL_FN_ { return v_[i]; }

    auto& weight() _PHARE_ALL_FN_ { return weight_; }
    auto& charge() _PHARE_ALL_FN_ { return charge_; }
    auto& iCell() _PHARE_ALL_FN_ { return iCell_; }
    auto& delta() _PHARE_ALL_FN_ { return delta_; }
    auto& v() _PHARE_ALL_FN_ { return v_; }

    auto& weight() const _PHARE_ALL_FN_ { return weight_; }
    auto& charge() const _PHARE_ALL_FN_ { return charge_; }
    auto& iCell() const _PHARE_ALL_FN_ { return iCell_; }
    auto& delta() const _PHARE_ALL_FN_ { return delta_; }
    auto& v() const _PHARE_ALL_FN_ { return v_; }


    // for performing the same operation across all vectors e.g. with std apply
    auto as_tuple(std::size_t i)
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple(std::size_t i) const
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }


    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args) _PHARE_ALL_FN_
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T*>{t, args...};
    }
    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, size()); }
    auto end() _PHARE_ALL_FN_ { return it(this, size()); }

    template<typename Impl, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void push_back(iterator_impl<Impl> const& particle)
    {
        this->weight_.push_back(particle.weight());
        this->charge_.push_back(particle.charge());
        this->iCell_.push_back(particle.iCell());
        this->delta_.push_back(particle.delta());
        this->v_.push_back(particle.v());
    }


    template<typename Particle_t, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void push_back(Particle_t const& particle)
    {
        auto const& [w, c, i, d, v, id] = particle;
        this->weight_.push_back(w);
        this->charge_.push_back(c);
        this->iCell_.push_back(i);
        this->delta_.push_back(d);
        this->v_.push_back(v);
    }

    template<typename Particle_t, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void emplace_back(Particle_t const& particle)
    {
        auto const& [w, c, i, d, v, id] = particle;
        this->weight_.push_back(w);
        this->charge_.push_back(c);
        this->iCell_.push_back(i);
        this->delta_.push_back(d);
        this->v_.push_back(v);
    }

    template<typename... Args>
    auto emplace_back(Args const&... args)
    {
        auto arg_tuple  = std::forward_as_tuple(args...);
        auto this_tuple = as_tuple();
        for_N<std::tuple_size_v<decltype(arg_tuple)>>([&](auto ic) {
            auto constexpr i = ic();
            std::get<i>(this_tuple).emplace_back(std::get<i>(arg_tuple));
        });

        return end() - 1;
    }


    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void clear()
    {
        std::apply([](auto&... container) { ((container.clear()), ...); }, as_tuple());
    }



    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void reserve(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.reserve(size)), ...); }, as_tuple());
    }

    void swap(This& that)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = that.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple).swap(std::get<i>(that_tuple)); });
    }

    void swap(std::size_t const& a, std::size_t const& b)
    {
        if (a == b)
            return;

        std::swap(weight_[a], weight_[b]);
        std::swap(charge_[a], charge_[b]);
        std::swap(iCell_[a], iCell_[b]);
        std::swap(delta_[a], delta_[b]);
        std::swap(v_[a], v_[b]);
    }


    // template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    // auto erase(iterator_impl<This*> a, iterator_impl<This*> b)
    // {
    //     return Super::erase(a, b);
    // }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void erase(iterator_impl<This*> first, iterator_impl<This*> last)
    {
        std::apply(
            [&](auto&... container) {
                ((container.erase(container.begin() + first.curr_pos,
                                  container.begin() + last.curr_pos)),
                 ...);
            },
            as_tuple());
    }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    std::size_t capacity() const
    {
        return weight_.capacity(); // they're all the same
    }


    auto copy(std::size_t i) const _PHARE_ALL_FN_
    {
        return Particle<dimension>{
            weight_[i], charge_[i], iCell_[i], delta_[i], v_[i],
        };
    }


    auto back() { return (*this)[size() - 1]; }
    auto front() { return (*this)[0]; }

    auto operator[](std::size_t const& s) const { return copy(s); }

    void replace_from(This const& that)
    {
        if (this == &that) // just in case
            return;
        auto src_tuple = that.as_tuple();
        auto dst_tuple = this->as_tuple();
        for_N<std::tuple_size_v<decltype(src_tuple)>>(
            [&](auto i) { std::get<i>(dst_tuple) = std::get<i>(src_tuple); });
    }

    auto constexpr static size_of_particle()
    {
        return sizeof(typename decltype(iCell_)::value_type)
               + sizeof(typename decltype(delta_)::value_type)
               + sizeof(typename decltype(weight_)::value_type)
               + sizeof(typename decltype(charge_)::value_type)
               + sizeof(typename decltype(v_)::value_type);
    }


    void check() const {}
};


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using SoAVectorParticles = SoAParticles<SoAVector<dim, alloc_mode>>;

template<std::size_t dim, std::size_t size>
using SoAArrayParticles = SoAParticles<SoAArray<dim, size>>;



template<typename OuterSuper>
template<typename T>
struct SoAParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;
    auto static constexpr is_const  = std::is_const_v<T>;
    auto static constexpr is_soa    = true;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = Particle<dimension>;
    using pointer           = Particle<dimension>*;
    using reference         = Particle<dimension>&;

    iterator_impl(T& particles_, std::size_t const& s = 0) _PHARE_ALL_FN_ : particles{particles_},
                                                                            curr_pos{s}
    {
    }
    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    auto& operator++()
    {
        ++curr_pos;
        return *this;
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto& operator+=(std::int64_t i) _PHARE_ALL_FN_
    {
        curr_pos += i;
        return *this;
    }


    auto& operator--()
    {
        --curr_pos;
        return *this;
    }
    auto operator+(std::int64_t i) const _PHARE_ALL_FN_
    {
        auto copy = *this;
        copy.curr_pos += i;
        return copy;
    }
    auto operator-(std::int64_t i) const
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
    auto operator-(iterator_impl const& that) { return curr_pos - that.curr_pos; }
    auto& operator()() { return deref(particles); }
    auto& operator()() const { return deref(particles); }
    auto& operator*() _PHARE_ALL_FN_ { return *this; }
    auto& operator*() const _PHARE_ALL_FN_ { return *this; }

    auto idx() const { return curr_pos; }

    auto& weight() _PHARE_ALL_FN_ { return deref(particles).weight_[curr_pos]; }
    auto& weight() const _PHARE_ALL_FN_ { return deref(particles).weight_[curr_pos]; }
    auto& charge() _PHARE_ALL_FN_ { return deref(particles).charge_[curr_pos]; }
    auto& charge() const _PHARE_ALL_FN_ { return deref(particles).charge_[curr_pos]; }
    auto& iCell() _PHARE_ALL_FN_ { return deref(particles).iCell_[curr_pos]; }
    auto& iCell() const _PHARE_ALL_FN_ { return deref(particles).iCell_[curr_pos]; }
    auto& delta() _PHARE_ALL_FN_ { return deref(particles).delta_[curr_pos]; }
    auto& delta() const _PHARE_ALL_FN_ { return deref(particles).delta_[curr_pos]; }
    auto& v() _PHARE_ALL_FN_ { return deref(particles).v_[curr_pos]; }
    auto& v() const _PHARE_ALL_FN_ { return deref(particles).v_[curr_pos]; }

    Particle<dimension> copy() const _PHARE_ALL_FN_
    {
        return {weight(), charge(), iCell(), delta(), v()};
    }

    T particles;
    std::size_t curr_pos = 0;
};


template<std::size_t dim, bool _const_ = false>
using ParticleArray_SOAView = SoAParticles<SoASpan<dim, _const_>>;


} // namespace PHARE::core


namespace std
{
template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_soa, PHARE::core::Particle<Iterator<O>::dimension>>
copy(Iterator<O> src)
{
    return {src.weight(), src.charge(), src.iCell(), src.delta(), src.v()};
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_soa, void> copy(Iterator<O> src_begin, //
                                                          Iterator<O> src_end,
                                                          Iterator<O> dst_begin)
{
    auto src_tuple = src_begin.particles.as_tuple();
    auto dst_tuple = dst_begin.particles.as_tuple();
    for_N<std::tuple_size_v<decltype(src_tuple)>>([&](auto i) {
        auto& src = std::get<i>(src_tuple);
        auto& dst = std::get<i>(dst_tuple);
        std::copy(&src[src_begin.curr_pos], &src[src_end.curr_pos], &dst[dst_begin.curr_pos]);
    });
}



template<template<typename> typename Iterator, typename O, typename Inserter>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_soa, void> copy(Iterator<O> src_begin, //
                                                          Iterator<O> src_end, Inserter i)
{
    // auto const& particles = src_begin.particles;

    for (; src_begin != src_end; ++src_begin, ++i)
        *i = src_begin.copy();
}

template<typename OuterSuper, typename T, typename Inserter, typename Fn,
         typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T>>
void copy_if(Iter src_begin, Iter src_end, Inserter i, Fn check)
{
    auto const& particles = src_begin.particles;

    for (; src_begin != src_end; ++src_begin)
        if (check(src_begin))
        {
            *i = src_begin.copy();
            ++i;
        }
}



template<typename OuterSuper, typename T,
         typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T>>
void iter_swap(Iter a, Iter b)
{
    PHARE_LOG_LINE_STR("");
    PHARE_ASSERT(&a.particles == &b.particles);

    a.particles.swap(a.curr_pos, b.curr_pos);
}

template<typename OuterSuper, typename T, typename Predicate,
         typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T>>
Iter find_if_not(Iter first, Iter last, Predicate q)
{ // https://en.cppreference.com/w/cpp/algorithm/find
    PHARE_LOG_LINE_STR("");
    for (; first != last; ++first)
        if (!q(first)) // pass iterator!
            return first;
    return last;
}

template<typename OuterSuper, typename T, typename Predicate/*,
         typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T>*/>
typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> partition(typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> first, typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> last, Predicate p)
{ // https://en.cppreference.com/w/cpp/algorithm/partition
    PHARE_LOG_LINE_STR("");
    first = find_if_not(first, last, p);
    if (first == last)
        return first;

    for (auto i = first + 1; i != last; ++i)
    {
        if (p(i)) // pass iterator!
        {
            iter_swap(i, first);
            ++first;
        }
    }
    return first;
}




template<template<typename> typename Iterator, typename O, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_soa, void>
transform(Iterator<O> src_begin, Iterator<O> src_end, Iterator<O> dst_begin, Fn transform)
{
    if (&src_begin.particles != &dst_begin.particles)
        throw std::runtime_error("transform cannot work");

    for (; src_begin != src_end; ++src_begin, ++dst_begin)
        /*dst_begin = */ transform(src_begin); // transform edits in place
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_soa, std::size_t> distance(Iterator<O> const& a,
                                                                     Iterator<O> const& b)
{
    PHARE_ASSERT(&deref(a.particles).weight(0) == &deref(b.particles).weight(0));
    PHARE_ASSERT(a.curr_pos <= b.curr_pos);

    return b.curr_pos - a.curr_pos;
}

} // namespace std


namespace PHARE::core
{


template<typename Super>
SoAParticles<Super>::SoAParticles(SoAParticles<Super>::iterator start,
                                  SoAParticles<Super>::iterator end)
    : Super{std::distance(start, end)}
{
    std::copy(start, end, this->begin()); // impl above
}

} // namespace PHARE::core


#if __has_include(<thrust/iterator/zip_iterator.h>)

#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/iterator/zip_iterator.h>

namespace PHARE::core::detail
{

struct SoAIteratorAdaptor
{
    template<typename SoAParticles_t>
    auto static make(SoAParticles_t& particles)
    {
        return thrust::make_zip_iterator( //
            thrust::make_tuple(           //
                particles.charge_.data(), // 0
                particles.weight_.data(), // 1
                particles.iCell_.data(),  // 2
                particles.delta_.data(),  // 3
                particles.v_.data()       // 4
                ));
    }

    template<typename T>
    static auto& iCell(T const& it) _PHARE_ALL_FN_
    {
        return thrust::get<2>(it);
    }
    template<typename T>
    static auto& delta(T const& it) _PHARE_ALL_FN_
    {
        return thrust::get<3>(it);
    }
};

} // namespace PHARE::core::detail


namespace PHARE::core
{


template<typename T, std::size_t dim>
auto partitionner(detail::SoAIteratorAdaptor& begin, detail::SoAIteratorAdaptor& end,
                  Box<T, dim> const& box)
{
    return partitionner(begin, end, box, [&box](auto const& part) {
        return isIn(detail::SoAIteratorAdaptor::iCell(part), box);
    });
}

} // namespace PHARE::core


#endif // if thrust header is found


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP */
