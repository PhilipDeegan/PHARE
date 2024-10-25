// C++11 version of https://en.cppreference.com/w/cpp/container/span

#ifndef PHARE_CORE_UTILITIES_SPAN_HPP
#define PHARE_CORE_UTILITIES_SPAN_HPP

#include <vector>
#include <cstddef>
#include <numeric>

#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/def.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

namespace PHARE::core
{
using span_size_default_t = unsigned long long int;

template<typename T, typename SIZE = span_size_default_t>
struct Span
{
    using value_type = std::decay_t<T>;

    Span(T* ptr_ = nullptr, SIZE s_ = 0)
        : ptr{ptr_}
        , s{s_}
    {
    }

    Span(Span&&)                 = default;
    Span(Span const&)            = default;
    Span& operator=(Span&&)      = default;
    Span& operator=(Span const&) = default;

    NO_DISCARD auto& operator[](SIZE i) _PHARE_ALL_FN_ { return ptr[i]; }
    NO_DISCARD auto& operator[](SIZE i) const _PHARE_ALL_FN_ { return ptr[i]; }
    NO_DISCARD T const* cdata() const _PHARE_ALL_FN_ { return ptr; }
    NO_DISCARD auto data() const _PHARE_ALL_FN_ { return ptr; }
    NO_DISCARD auto data() _PHARE_ALL_FN_ { return ptr; }
    // NO_DISCARD T const* const& data() const { return ptr; }
    NO_DISCARD auto begin() _PHARE_ALL_FN_ { return ptr; }
    NO_DISCARD auto begin() const _PHARE_ALL_FN_ { return ptr; }
    NO_DISCARD auto end() _PHARE_ALL_FN_ { return ptr + s; }
    NO_DISCARD auto end() const _PHARE_ALL_FN_ { return ptr + s; }
    NO_DISCARD SIZE const& size() const _PHARE_ALL_FN_ { return s; }

    T* ptr = nullptr;
    SIZE s = 0;
};

template<typename T, typename data_fn = void, typename size_fn = void>
struct is_span_like : std::false_type
{
    // not a span
};

template<typename T>
struct is_span_like<T, core::tryToInstanciate<decltype(std::declval<T>().data())>,
                    core::tryToInstanciate<decltype(std::declval<T>().size())>> : std::true_type
{
    // is interopable with span i.e. has a data() and size() function
};

template<typename T>
auto constexpr is_span_like_v = is_span_like<T>::value;

template<typename Container, std::enable_if_t<core::is_span_like_v<Container>, bool> = 0>
auto make_span(Container& container)
{
    return Span<typename Container::value_type>{container.data(), container.size()};
}
template<typename Container, std::enable_if_t<core::is_span_like_v<Container>, bool> = 0>
auto make_span(Container& container, std::size_t const& size)
{
    return Span<typename Container::value_type>{container.data(), size};
}
template<typename Container, std::enable_if_t<core::is_span_like_v<Container>, bool> = 0>
auto make_span(Container const& container)
{
    return Span<typename Container::value_type>{container.data(), container.size()};
}


template<typename T, typename SIZE = span_size_default_t>
class VectorSpan : private StackVar<std::vector<T>>, public core::Span<T, SIZE>
{
    using Vector = StackVar<std::vector<T>>;
    using Span_  = Span<T, SIZE>;

public:
    VectorSpan(SIZE size, T value)
        : Vector{std::vector<T>(size, value)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
    VectorSpan(std::vector<T>&& vec_)
        : Vector{std::move(vec_)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
    VectorSpan(std::vector<T> const& vec_)
        : Vector{vec_}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
};



template<typename T, typename SIZE = span_size_default_t>
struct SpanSet
{
    using value_type = T;
    using SpanSet_   = SpanSet<T, SIZE>;

    SpanSet() = default;

    SpanSet(std::vector<SIZE>&& sizes_)
        : size{std::accumulate(sizes_.begin(), sizes_.end(), 0)}
        , sizes(sizes_)
        , displs(core::displacementFrom(sizes))
        , vec(size)
    {
    }

    SpanSet(SpanSet&& from)
        : size{from.size}
        , sizes(std::move(from.sizes))
        , displs(std::move(from.displs))
        , vec(std::move(from.vec))
    {
    }


    NO_DISCARD Span<T, SIZE> operator[](SIZE i)
    {
        return {this->vec.data() + displs[i], this->sizes[i]};
    }
    NO_DISCARD Span<T const, SIZE> operator[](SIZE i) const
    {
        return {this->vec.data() + displs[i], this->sizes[i]};
    }

    NO_DISCARD auto data() const { return vec.data(); }
    NO_DISCARD auto data() { return const_cast<T*>(vec.data()); }

    struct iterator
    {
        iterator(SpanSet_* _sv)
            : sv(_sv)
        {
        }
        iterator operator++()
        {
            curr_pos += sv->sizes[curr_ptr++];
            return *this;
        }
        bool operator!=(iterator const& /*other*/) const { return curr_ptr != sv->sizes.size(); }
        Span<T, SIZE> operator*() const { return {sv->vec.data() + curr_pos, sv->sizes[curr_ptr]}; }

        SpanSet_* sv  = nullptr;
        SIZE curr_pos = 0, curr_ptr = 0;
    };

    NO_DISCARD auto begin() { return iterator(this); }
    NO_DISCARD auto cbegin() const { return iterator(this); }

    NO_DISCARD auto end() { return iterator(this); }
    NO_DISCARD auto cend() const { return iterator(this); }

    SIZE size;
    std::vector<SIZE> sizes;
    std::vector<SIZE> displs;
    std::vector<T> vec;
};




template<typename T, std::size_t size>
auto flatten(std::vector<std::array<T, size>> const& data)
{
    assert(data.size() > 0);

    return Span<T const, std::size_t>{data.data()->data(), data.size() * size};
}


template<typename T, std::size_t size>
auto flatten(std::vector<std::array<T, size>>& data)
{
    assert(data.size() > 0);

    return Span<T, std::size_t>{data.data()->data(), data.size() * size};
}



} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_SPAN_HPP
