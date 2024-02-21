#ifndef PHARE_PYTHON_INTEPRETER_HPP
#define PHARE_PYTHON_INTEPRETER_HPP

#include "core/logger.hpp"

// pybind or not
#if __has_include("pybind11/embed.h")

// clang-format off
DISABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)

#undef HAVE_SYS_TIMES_H // included in python again, possibly with different value
#undef HAVE_UNISTD_H

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/functional.h>
ENABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)
// clang-format on

#else // assume nanobind for now

// #include "numpy/arrayobject.h"
// #include "numpy/npy_math.h"



int _init_python()
{
    PyStatus status;
    PyConfig config;
    PyConfig_InitPythonConfig(&config);
    status = PyConfig_SetString(&config, &config.program_name, L"/path/to/my_program");
    if (PyStatus_Exception(status))
        throw std::runtime_error("no");
    status = Py_InitializeFromConfig(&config);
    if (PyStatus_Exception(status))
        throw std::runtime_error("no");
    return 0;
}
static const int init_python_ = _init_python();

#endif




namespace PHARE::py3
{


class Interpreter
{
public:
    static Interpreter& INSTANCE()
    {
        static Interpreter i;
        return i;
    }

    Interpreter() {}
    ~Interpreter() { Py_Finalize(); }

    void import()
    {
#if __has_include("pybind11/embed.h")
        auto module = py::module::import(initModuleName_.c_str());
        module.attr(functionName_)(moduleName_);
#else

        auto module    = nanobind::module_::import_(initModuleName_.c_str());
        PyObject* name = PyUnicode_FromString(functionName_.c_str());
        PyObject* arg  = PyUnicode_FromString(moduleName_.c_str());
        auto res       = PyObject_CallMethodOneArg(module.ptr(), name, arg);
#endif
    }




private:
    std::string functionName_{"get_user_inputs"};
    std::string moduleName_{"job"};
    std::string initModuleName_{"pyphare.pharein.init"};
#if __has_include("pybind11/embed.h")
    py::scoped_interpreter guard_;
#else

#endif
};




} // namespace PHARE::py3

#endif /*PHARE_PYTHON_INTEPRETER_HPP*/
