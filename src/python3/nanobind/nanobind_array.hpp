#ifndef PHARE_PYTHON_NANOBIND_NANOBIND_ARRAY_HPP
#define PHARE_PYTHON_NANOBIND_NANOBIND_ARRAY_HPP

// interop with pybind mostly


#include "nanobind/ndarray.h"



namespace PHARE::pydata
{



template<typename T>
struct __attribute__((visibility("hidden"))) py_array_t
{
    py_array_t(std::size_t size)
        : ptr{reinterpret_cast<T*>(std::malloc(size * sizeof(T)))}
        , deleter{ptr, [](void* p) noexcept { std::free(p); }}
        , array{ptr, {size}, deleter}
    {
    }

    auto& request() { return *this; }

    T* ptr = nullptr;
    nanobind::capsule deleter;
    nanobind::ndarray<T, nanobind::c_contig> array;
};

// template<typename T>
// using py_array_t = nanobind::ndarray<T, nanobind::c_contig>;




} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_NANOBIND_NANOBIND_ARRAY_HPP*/
