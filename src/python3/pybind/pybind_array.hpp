#ifndef PHARE_PYTHON_NANOBIND_PYBIND_ARRAY_HPP
#define PHARE_PYTHON_NANOBIND_PYBIND_ARRAY_HPP


#include "pybind11/stl.h"
#include "pybind11/numpy.h"


namespace PHARE::pydata
{

template<typename T>
using py_array_t = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_NANOBIND_PYBIND_ARRAY_HPP*/
