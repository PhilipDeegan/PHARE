#ifndef PHARE_PYTHON3_TAGGING_HPP
#define PHARE_PYTHON3_TAGGING_HPP

// cpp_etc's tagging outside a simulation: a run's tagging strategy on fields given as numpy
// arrays (pyphare.simulator.tagger.Tagger)

#include "simulator/simulator_def.hpp"
#include "initializer/data_provider.hpp"
#include "amr/tagging/tagging_fields_model.hpp"

#include "python3/pybind_def.hpp"

#include <span>
#include <array>
#include <string>
#include <cstddef>
#include <stdexcept>

namespace PHARE::pydata
{
namespace py = pybind11;

/// a python dict of dicts, bools, ints, floats and strings as a PHAREDict
inline void toPHAREDict(py::dict const& in, initializer::PHAREDict& out)
{
    for (auto const& [k, v] : in)
    {
        auto const key = py::cast<std::string>(k);
        if (py::isinstance<py::dict>(v))
            toPHAREDict(py::cast<py::dict>(v), out[key]);
        else if (py::isinstance<py::bool_>(v))
            out[key] = py::cast<bool>(v);
        else if (py::isinstance<py::int_>(v))
            out[key] = py::cast<int>(v);
        else if (py::isinstance<py::float_>(v))
            out[key] = py::cast<double>(v);
        else if (py::isinstance<py::str>(v))
            out[key] = py::cast<std::string>(v);
        else
            throw std::runtime_error("tagging: dict value of '" + key
                                     + "' is not a dict, bool, int, float or string");
    }
}

template<std::size_t dim, std::size_t interp>
void declareTagPatch(py::module& m)
{
    auto constexpr opts    = SimOpts{.dimension = dim, .interp_order = interp};
    std::string const name = "tag_patch_" + std::to_string(dim) + "_" + std::to_string(interp);

    m.def(
        name.c_str(),
        [](py::dict const& tagging, std::array<double, dim> const& meshSize,
           std::array<int, dim> const& lower, std::array<int, dim> const& upper,
           py::array_t<double, py::array::c_style> bx, py::array_t<double, py::array::c_style> by,
           py::array_t<double, py::array::c_style> bz) {
            initializer::PHAREDict dict;
            toPHAREDict(tagging, dict);
            auto span = [](auto& a) {
                return std::span<double>{a.mutable_data(), static_cast<std::size_t>(a.size())};
            };
            auto const tags = amr::tagging::tagPatch<opts>(dict, meshSize, lower, upper,
                                                           {span(bx), span(by), span(bz)});
            return py::array_t<int>(static_cast<py::ssize_t>(tags.size()), tags.data());
        },
        "the tags (Fortran ordering) of the patch lower..upper from B over its allocation, with "
        "the tagging config dict a run reads at simulation/AMR/refinement/tagging");
}

inline void declareTagging(py::module& m)
{
    declareTagPatch<1, 1>(m);
    declareTagPatch<1, 2>(m);
    declareTagPatch<1, 3>(m);
    declareTagPatch<2, 1>(m);
    declareTagPatch<2, 2>(m);
    declareTagPatch<2, 3>(m);
    declareTagPatch<3, 1>(m);
    declareTagPatch<3, 2>(m);
    declareTagPatch<3, 3>(m);
}

} // namespace PHARE::pydata

#endif
