
// Philip Deegan comments: IWYU may work on this file, but it takes forever so I never found out.
//   headers were included manually

#include <vector>
#include <cstddef>

#include "mpi.h"

#include "core/utilities/mpi_utils.h"
#include "core/data/particles/particle.h"
#include "core/data/particles/contiguous.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/types.h"
#include "core/utilities/mpi_utils.h"
#include "amr/wrappers/hierarchy.h"
#include "phare/phare.h"
#include "simulator/simulator.h"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include "python3/pybind_def.h"
#include "python3/particles.h"
#include "python3/patch_data.h"
#include "python3/patch_level.h"
#include "python3/data_wrangler.h"


#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME cpp
#endif


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

template<typename Float>
std::string float_id()
{
    if constexpr (std::is_same_v<Float, double>)
        return "d";
    else if constexpr (std::is_same_v<Float, float>)
        return "s";
    else
        throw std::runtime_error("Unhandled data type");
}

template<typename Float, std::size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<Float, dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim) + "_" + float_id<Float>();
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

template<typename _dim, typename _interp, typename _nbRefinedPart, typename Float>
void declare_sim(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();


    std::string type_str = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                           + std::to_string(nbRefinedPart);
    std::string type_string = type_str + "_" + float_id<Float>();

    using Sim        = Simulator<dim, interp, nbRefinedPart, Float>;
    std::string name = "Simulator" + type_string;
    declareSimulator<Sim>(
        py::class_<Sim, std::shared_ptr<Sim>>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));


    name = "make_simulator" + type_string;
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy<Float>> const& hier) {
        return std::shared_ptr<Sim>{
            std::move(makeSimulator<dim, interp, nbRefinedPart, Float>(hier))};
    });


    using DW = DataWrangler<dim, interp, nbRefinedPart, Float>;
    name     = "DataWrangler" + type_string;
    py::class_<DW, std::shared_ptr<DW>>(m, name.c_str())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy<Float>> const&>())
        .def(py::init<std::shared_ptr<ISimulator<Float>> const&,
                      std::shared_ptr<amr::Hierarchy<Float>> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getPatchLevel", &DW::getPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    using PL = PatchLevel<dim, interp, nbRefinedPart, Float>;
    name     = "PatchLevel_" + type_string;

    py::class_<PL, std::shared_ptr<PL>>(m, name.c_str())
        .def("getEM", &PL::getEM)
        .def("getE", &PL::getE)
        .def("getB", &PL::getB)
        .def("getBx", &PL::getBx)
        .def("getBy", &PL::getBy)
        .def("getBz", &PL::getBz)
        .def("getEx", &PL::getEx)
        .def("getEy", &PL::getEy)
        .def("getEz", &PL::getEz)
        .def("getVix", &PL::getVix)
        .def("getViy", &PL::getViy)
        .def("getViz", &PL::getViz)
        .def("getDensity", &PL::getDensity)
        .def("getBulkVelocity", &PL::getBulkVelocity)
        .def("getPopDensities", &PL::getPopDensities)
        .def("getPopFluxes", &PL::getPopFlux)
        .def("getFx", &PL::getFx)
        .def("getFy", &PL::getFy)
        .def("getFz", &PL::getFz)
        .def("getParticles", &PL::getParticles, py::arg("userPopName") = "all");

    using _Splitter
        = PHARE::amr::Splitter<_dim, _interp, core::RefinedParticlesConst<nbRefinedPart>>;
    if constexpr (std::is_same_v<Float, double>) // only needs one declarion for both double/float
    {
        name = "Splitter" + type_str;
        py::class_<_Splitter, std::shared_ptr<_Splitter>>(m, name.c_str())
            .def(py::init<>())
            .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
            .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });
    }

    name = "split_pyarray_particles" + type_string;
    m.def(name.c_str(), splitPyArrayParticles<_Splitter, Float>);
}


template<typename Dimension, typename InterpOrder, typename Float, typename... NbRefinedParts>
void declare(py::module& m, std::tuple<Dimension, InterpOrder, Float, NbRefinedParts...> const&)
{
    core::apply(std::tuple<NbRefinedParts...>{}, [&](auto& nbRefinedPart) {
        declare_sim<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>, Float>(m);
    });
}



auto pybind_version()
{
    std::stringstream ss;
    ss << PYBIND11_VERSION_MAJOR << ".";
    ss << PYBIND11_VERSION_MINOR << ".";
    ss << PHARE_TO_STR(PYBIND11_VERSION_PATCH);
    return ss.str();
}

auto samrai_version()
{
    std::stringstream ss;
    ss << SAMRAI_VERSION_MAJOR << ".";
    ss << SAMRAI_VERSION_MINOR << ".";
    ss << SAMRAI_VERSION_PATCHLEVEL;
    return ss.str();
}



PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m)
{
    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy<double>, std::shared_ptr<PHARE::amr::Hierarchy<double>>>(
        m, "AMRHierarchy_d");
    // py::class_<PHARE::amr::Hierarchy<float>, std::shared_ptr<PHARE::amr::Hierarchy<float>>>(
    //     m, "AMRHierarchy_s");

    declareSimulator<ISimulator<double>>(
        py::class_<ISimulator<double>, std::shared_ptr<ISimulator<double>>>(m, "ISimulator_d")
            .def("interp_order", &ISimulator<double>::interporder)
            .def("dump", &ISimulator<double>::dump, py::arg("timestamp"), py::arg("timestep")));

    // declareSimulator<ISimulator<float>>(
    //     py::class_<ISimulator<float>, std::shared_ptr<ISimulator<float>>>(m, "ISimulator_s")
    //         .def("interp_order", &ISimulator<float>::interporder)
    //         .def("dump", &ISimulator<float>::dump, py::arg("timestamp"), py::arg("timestep")));

    m.def("make_hierarchy_d", []() { return PHARE::amr::Hierarchy<double>::make(); });
    // m.def("make_hierarchy_s", []() { return PHARE::amr::Hierarchy<float>::make(); });

    m.def("mpi_size", []() { return core::mpi::size(); });
    m.def("mpi_rank", []() { return core::mpi::rank(); });

    declareDim<double, 1>(m);
    declareDim<double, 2>(m);
    declareDim<double, 3>(m);

    // declareDim<float, 1>(m);
    // declareDim<float, 2>(m);
    // declareDim<float, 3>(m);

    core::apply(core::possibleSimulators<double>(),
                [&](auto const& simType) { declare(m, simType); });
    // core::apply(core::possibleSimulators<float>(),
    //             [&](auto const& simType) { declare(m, simType); });

    declarePatchData<std::vector<double>, 1>(m, "PatchDataVectorDouble_1D");
    declarePatchData<std::vector<double>, 2>(m, "PatchDataVectorDouble_2D");
    declarePatchData<std::vector<double>, 3>(m, "PatchDataVectorDouble_3D");

    py::class_<core::Span<double>, std::shared_ptr<core::Span<double>>>(m, "Span_d");
    // py::class_<core::Span<float>, std::shared_ptr<core::Span<float>>>(m, "Span_s");

    py::class_<PyArrayWrapper<double>, std::shared_ptr<PyArrayWrapper<double>>, core::Span<double>>(
        m, "PyWrapper_d");
    // py::class_<PyArrayWrapper<float>, std::shared_ptr<PyArrayWrapper<float>>, core::Span<float>>(
    //     m, "PyWrapper_s");

    m.def("phare_deps", []() {
        std::unordered_map<std::string, std::string> versions{{"pybind", pybind_version()},
                                                              {"samrai", samrai_version()}};
        _PHARE_WITH_HIGHFIVE(versions["highfive"] = PHARE_TO_STR(HIGHFIVE_VERSION));
        return versions;
    });

    m.def("makePyArrayWrapper_d", makePyArrayWrapper<double>);
    // m.def("makePyArrayWrapper_s", makePyArrayWrapper<float>);
}


} // namespace PHARE::pydata
