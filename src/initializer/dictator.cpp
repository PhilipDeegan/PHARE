
#include "cppdict/include/dict.hpp"
#include "initializer/python_data_provider.h"

#include "python3/pybind_def.h"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>



namespace py = pybind11;

using PHARE::initializer::InitFunction;



template<typename T>
void add(std::string const& path, T&& value)
{
    cppdict::add(path, std::forward<T>(value),
                 PHARE::initializer::PHAREDictHandler::INSTANCE().dict());
}

template<typename T>
void add_array_as_vector(std::string path, PHARE::pydata::py_array_t<T>& array)
{
    auto buf = array.request();

    if (buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    T* ptr = static_cast<T*>(buf.ptr);

    cppdict::add(path, std::vector<T>(ptr, ptr + buf.size),
                 PHARE::initializer::PHAREDictHandler::INSTANCE().dict());
}


PYBIND11_MODULE(dictator, m)
{
    // expose dict add function per template type to force casts/error when used from python
    m.def("add_size_t", add<std::size_t>, "add_size_t");
    m.def("add_optional_size_t", add<std::optional<std::size_t>>, "add_optional_size_t");

    m.def("add_int", add<int>, "add");
    m.def("add_vector_int", add<std::vector<int>>, "add");
    m.def("add_float", add<float>, "add");
    m.def("add_double", add<double>, "add");
    m.def("add_s", add<float>, "add");
    m.def("add_d", add<double>, "add");
    m.def("add_string", add<std::string>, "add");

    m.def("addInitFunction_d_1D", add<InitFunction<double, 1>>, "add");
    m.def("addInitFunction_d_2D", add<InitFunction<double, 2>>, "add");
    m.def("addInitFunction_d_3D", add<InitFunction<double, 3>>, "add");

    m.def("addInitFunction_s_1D", add<InitFunction<float, 1>>, "add");
    m.def("addInitFunction_s_2D", add<InitFunction<float, 2>>, "add");
    m.def("addInitFunction_s_3D", add<InitFunction<float, 3>>, "add");

    m.def("add_array_as_vector", add_array_as_vector<double>, "add_array_as_vector");
}
