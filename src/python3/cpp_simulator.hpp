#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP

#include <vector>
#include <cstddef>

#include "core/def/phare_mpi.hpp"

#include "core/utilities/mpi_utils.hpp"
#include "core/data/particles/particle.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "phare/phare.hpp"
#include "simulator/simulator.hpp"

#include "python3/pybind_def.hpp"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include "python3/particles.hpp"
#include "python3/patch_data.hpp"
#include "python3/patch_level.hpp"
#include "python3/data_wrangler.hpp"



namespace py = pybind11;

namespace PHARE::pydata
{
template<typename Type, std::size_t dimension>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type, dimension>;
    py::class_<PatchDataType>(m, key.c_str())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("lower", &PatchDataType::lower)
        .def_readonly("upper", &PatchDataType::upper)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<std::size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim);
    py::class_<CP, std::shared_ptr<CP>>(m, name.c_str())
        .def(py::init<std::size_t>())
        .def_readwrite("iCell", &CP::iCell)
        .def_readwrite("delta", &CP::delta)
        .def_readwrite("weight", &CP::weight)
        .def_readwrite("charge", &CP::charge)
        .def_readwrite("v", &CP::v)
        .def("size", &CP::size);

    name = "PatchData" + name;
    declarePatchData<CP, dim>(m, name.c_str());
}

template<typename Simulator, typename PyClass>
void declareSimulator(PyClass&& sim)
{
    sim.def("initialize", &Simulator::initialize)
        .def("advance", &Simulator::advance)
        .def("startTime", &Simulator::startTime)
        .def("currentTime", &Simulator::currentTime)
        .def("endTime", &Simulator::endTime)
        .def("timeStep", &Simulator::timeStep)
        .def("to_str", &Simulator::to_str)
        .def("domain_box", &Simulator::domainBox)
        .def("cell_width", &Simulator::cellWidth)
        .def("dump", &Simulator::dump, py::arg("timestamp"), py::arg("timestep"));
}

template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare_etc(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();

    std::string const dimterp_string = "_" + std::to_string(dim) + "_" + std::to_string(interp);
    std::string const type_string    = dimterp_string + "_" + std::to_string(nbRefinedPart);

    using Sim        = Simulator<dim, interp, nbRefinedPart>;
    using DW         = DataWrangler<dim, interp, nbRefinedPart>;
    std::string name = "DataWrangler" + type_string;
    py::class_<DW, std::shared_ptr<DW>>(m, name.c_str())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def(py::init<std::shared_ptr<ISimulator> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getMHDPatchLevel", &DW::getMHDPatchLevel)
        .def("getHybridPatchLevel", &DW::getHybridPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    using HybPL = PatchLevel<typename Sim::HybridModel>;
    name        = "HybridPatchLevel_" + dimterp_string;
    if constexpr (core::defaultNbrRefinedParts<dim, interp>() == nbRefinedPart) // register once!
        py::class_<HybPL, std::shared_ptr<HybPL>>(m, name.c_str())
            .def("getEM", &HybPL::getEM)
            .def("getE", &HybPL::getE)
            .def("getB", &HybPL::getB)
            .def("getBx", &HybPL::getBx)
            .def("getBy", &HybPL::getBy)
            .def("getBz", &HybPL::getBz)
            .def("getEx", &HybPL::getEx)
            .def("getEy", &HybPL::getEy)
            .def("getEz", &HybPL::getEz)
            .def("getVi", &HybPL::getVi)
            .def("getVix", &HybPL::getVix)
            .def("getViy", &HybPL::getViy)
            .def("getViz", &HybPL::getViz)
            .def("getDensity", &HybPL::getDensity)
            .def("getBulkVelocity", &HybPL::getBulkVelocity)
            .def("getPopDensities", &HybPL::getPopDensities)
            .def("getPopFluxes", &HybPL::getPopFlux)
            .def("getFx", &HybPL::getFx)
            .def("getFy", &HybPL::getFy)
            .def("getFz", &HybPL::getFz)
            .def("getParticles", &HybPL::getParticles, py::arg("userPopName") = "all");

    using MHDPL = PatchLevel<typename Sim::MHDModel>;
    name        = "MHDPatchLevel_" + dimterp_string;
    if constexpr (core::defaultNbrRefinedParts<dim, interp>() == nbRefinedPart)
        py::class_<MHDPL, std::shared_ptr<MHDPL>>(m, name.c_str());

    using _Splitter
        = PHARE::amr::Splitter<_dim, _interp, core::RefinedParticlesConst<nbRefinedPart>>;
    name = "Splitter" + type_string;
    py::class_<_Splitter, std::shared_ptr<_Splitter>>(m, name.c_str())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

    name = "split_pyarray_particles" + type_string;
    m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
}

template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare_sim(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();

    std::string type_string = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                              + std::to_string(nbRefinedPart);

    using Sim        = Simulator<dim, interp, nbRefinedPart>;
    std::string name = "Simulator" + type_string;
    declareSimulator<Sim>(
        py::class_<Sim, std::shared_ptr<Sim>>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));

    name = "make_simulator" + type_string;
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::shared_ptr<Sim>{std::move(makeSimulator<dim, interp, nbRefinedPart>(hier))};
    });
}

template<typename Dimension, typename InterpOrder, typename... NbRefinedParts>
void declare_all(py::module& m, std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
{
    core::apply(std::tuple<NbRefinedParts...>{}, [&](auto& nbRefinedPart) {
        declare_sim<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
        declare_etc<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
    });
}

void inline declare_essential(py::module& m)
{
    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy")
        .def("modelForLevel", &PHARE::amr::Hierarchy::modelForLevel);

    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });

    m.def("mpi_size", []() { return core::mpi::size(); });
    m.def("mpi_rank", []() { return core::mpi::rank(); });
    m.def("mpi_barrier", []() { core::mpi::barrier(); });
}


} // namespace PHARE::pydata


// https://stackoverflow.com/a/51061314/795574
// ASAN detects leaks by default, even in system/third party libraries
inline char const* __asan_default_options()
{
    return "detect_leaks=0";
}



#endif /*PHARE_PYTHON_CPP_SIMULATOR_H*/
