#ifndef PHARE_PYTHON_PYBIND_DEF_H
#define PHARE_PYTHON_PYBIND_DEF_H

#include <tuple>
#include <cstdint>
#include <stdexcept>

#include "core/utilities/span.h"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"

namespace PHARE::pydata
{
template<typename T>
using py_array_t = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


using pyarray_particles_t = std::tuple<py_array_t<int32_t>, py_array_t<float>, py_array_t<double>,
                                       py_array_t<double>, py_array_t<double>>;

using pyarray_particles_crt
    = std::tuple<py_array_t<int32_t> const&, py_array_t<float> const&, py_array_t<double> const&,
                 py_array_t<double> const&, py_array_t<double> const&>;

template<typename PyArrayInfo>
std::size_t ndSize(PyArrayInfo const& ar_info)
{
    assert(ar_info.ndim >= 1 && ar_info.ndim <= 3);

    std::size_t size = ar_info.shape[0];
    for (size_t i = 1; i < ar_info.ndim; i++)
        size *= ar_info.shape[i];

    return size;
}


template<typename T>
class PyArrayWrapper : public core::Span<T>
{
public:
    PyArrayWrapper(PHARE::pydata::py_array_t<T> const& array)
        : core::Span<T>{static_cast<T*>(array.request().ptr), pydata::ndSize(array.request())}
        , _array{array}
    {
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
    return {static_cast<T*>(ar_info.ptr), ndSize(ar_info)};
}



template<typename T>
core::Span<T, int> to_span(py_array_t<T> const& py_array)
{
    py::buffer_info info = py_array.request();
    if (!info.ptr)
        throw std::runtime_error("to_span received an array with an invalid internal ptr");
    assert(info.ndim == 1 or info.ndim == 2);
    int size = info.ndim == 1 ? info.shape[0] : info.shape[0] * info.shape[1];
    return {static_cast<T*>(info.ptr), size};
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PYBIND_DEF_H*/
