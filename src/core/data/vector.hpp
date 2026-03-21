#ifndef PHARE_CORE_DATA_VECTOR_HPP
#define PHARE_CORE_DATA_VECTOR_HPP

// #include "core/def/phare_config.hpp"
// #include "core/utilities/types.hpp"

#include "core/utilities/allocators.hpp"

#include <vector>
#include <cstdint>

namespace PHARE::core
{


template<typename T, typename Allocator = NonConstructingAllocator<T>>
struct MinimizingVector
{
    using vector_t = std::vector<T, Allocator>;

    template<bool copy_old = true>
    auto& get(std::size_t const& s)
    {
        if (s < vec.capacity() * percentile)
            ++_c;
        else
            _c = 0;

        if (_c == period)
        {
            vector_t r;
            r.reserve(vec.capacity() * realloc_to);
            if constexpr (copy_old)
                r = vec;
            vec = std::move(r);
            _c  = 0;
        }

        vec.resize(s);
        return *this;
    }

    auto& get_no_copy(std::size_t const s) { return get<false>(s); }
    auto& resize(std::size_t const& s) { return get(s); }
    auto& clear()
    {
        vec.clear();
        return *this;
    }
    auto& reserve_and_clear(std::size_t const s) { return get<false>(s).clear(); }

    auto& operator[](std::size_t const i) { return vec.data()[i]; }
    auto& operator[](std::size_t const i) const { return vec.data()[i]; }
    auto& operator*() { return vec; }
    auto& operator()() { return vec; }
    auto& operator()() const { return vec; }

    auto& emplace_back(auto&&... args) { return vec.emplace_back(args...); }
    void push_back(auto&&... args) { vec.push_back(args...); }
    auto data() { return vec.data(); }
    auto data() const { return vec.data(); }
    auto size() const { return vec.size(); }
    auto capacity() const { return vec.capacity(); }

    void destroy() { vec = std::move(vector_t{Allocator{}}); }

    double const percentile  = .80;
    double const realloc_to  = .85;
    std::size_t const period = 100;
    vector_t vec{Allocator{}};
    std::uint16_t _c = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_VECTOR_HPP */