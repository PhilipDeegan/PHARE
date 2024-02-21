#ifndef PHARE_PYTHON_INTEROP_HPP
#define PHARE_PYTHON_INTEROP_HPP

#include <tuple>
#include <cassert>
#include <cstdint>
#include <stdexcept>

#include "core/utilities/span.hpp"

#undef HAVE_SYS_TIMES_H // included in python again, possibly with different value
#undef HAVE_UNISTD_H

#if __has_include("pybind11/stl.h")
#include "python3/pybind/pybind_array.hpp"

#else // assume nanobind for now
#include "python3/nanobind/nanobind_array.hpp"

#endif


namespace PHARE::pydata
{


using pyarray_particles_t = std::tuple<py_array_t<std::int32_t>, py_array_t<double>,
                                       py_array_t<double>, py_array_t<double>, py_array_t<double>>;

using pyarray_particles_crt
    = std::tuple<py_array_t<std::int32_t> const&, py_array_t<double> const&,
                 py_array_t<double> const&, py_array_t<double> const&, py_array_t<double> const&>;

template<typename PyArrayInfo>
std::size_t ndSize(PyArrayInfo const& ar_info)
{
    assert(ar_info.ndim >= 1 && ar_info.ndim <= 3);

    return std::accumulate(ar_info.shape.begin(), ar_info.shape.end(), 1,
                           std::multiplies<std::size_t>());
}


template<typename T>
class PyArrayWrapper : public core::Span<T>
{
public:
    PyArrayWrapper(PHARE::pydata::py_array_t<T> const& array)
        : core::Span<T>{static_cast<T*>(array.request().ptr), pydata::ndSize(array.request())}
        , _array{array}
    {
        assert(_array.request().ptr);
        assert(_array.request().ptr == array.request().ptr); // assert no copy
    }

protected:
    PHARE::pydata::py_array_t<T> _array;
};

template<typename T>
std::shared_ptr<core::Span<T>> makePyArrayWrapper(py_array_t<T> const& array)
{
    return std::make_shared<PyArrayWrapper<T>>(array);
}

template<typename T>
core::Span<T> makeSpan(py_array_t<T> const& py_array)
{
    auto ar_info = py_array.request();
    assert(ar_info.ptr);
    return {static_cast<T*>(ar_info.ptr), ndSize(ar_info)};
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PYBIND_DEF_H*/
