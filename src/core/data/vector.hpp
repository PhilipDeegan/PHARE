#ifndef PHARE_CORE_DATA_VECTOR_HPP
#define PHARE_CORE_DATA_VECTOR_HPP

// #include "core/def/phare_config.hpp"
// #include "core/utilities/types.hpp"

#include <vector>
#include <cstdint>

namespace PHARE::core
{


template<typename T>
struct MinimizingVector
{
    using vector_t = std::vector<T>;

    MinimizingVector() = default;
    MinimizingVector(std::size_t const s)
        : v(s)
    {
    }

    template<bool copy_old = true>
    auto& get(std::size_t const& s)
    {
        if (s < v.capacity() * percentile)
            ++_c;
        else
            _c = 0;

        if (_c == period)
        {
            vector_t r;
            r.reserve(v.capacity() * realloc_to);
            if constexpr (copy_old)
                r = v;
            v  = std::move(r);
            _c = 0;
        }

        v.resize(s);
        return v;
    }

    auto& get_no_copy(std::size_t const s) { return get<false>(s); }
    auto& resize(std::size_t const& s) { return get(s); }
    auto& reserve(std::size_t const& s)
    {
        get(s).resize(0);
        return v;
    }

    auto& zero(std::size_t const& s)
    {
        clear();

        if (s < size())
        {
            if (s < v.capacity() * percentile)
                ++_c;
            else
                _c = 0;

            if (_c == period)
            {
                vector_t r;
                r.reserve(v.capacity() * realloc_to);
                v  = std::move(r);
                _c = 0;
            }
        }
        else
            v.reserve(s);

        return v;
    }


    auto& operator[](std::size_t const i) { return v.data()[i]; }
    auto& operator[](std::size_t const i) const { return v.data()[i]; }
    auto& operator()() { return v; }
    auto& operator()() const { return v; }

    auto& clear()
    {
        if (v.size())
            v.clear();
        return v;
    }
    void push_back(auto&&... args) { v.push_back(args...); }
    auto size() const { return v.size(); }
    auto capacity() const { return v.capacity(); }

    double const percentile  = .80;
    double const realloc_to  = .85;
    std::size_t const period = 100;
    vector_t v;
    std::uint16_t _c = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_VECTOR_HPP */