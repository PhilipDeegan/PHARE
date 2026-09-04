// This file is for the python module for everything besides C++ Simulators.


#include "core/def.hpp"
#include "core/def/phare_config.hpp" // IWYU pragma: keep

#include "amr/samrai.hpp"             // SamraiLifeCycle without simulators
#include "amr/wrappers/hierarchy.hpp" // for HierarchyRestarter::getRestartFileFullPath

#include "python3/pybind_def.hpp"

#include "hdf5/phare_hdf5.hpp"

#include "simulator/simulator_runtime.hpp"

#if PHARE_HAS_HIGHFIVE
#include "hdf5/detail/h5/h5_file.hpp"
#endif

#include <pybind11/stl_bind.h>
#include <pybind11/native_enum.h>

namespace py = pybind11;

namespace PHARE::pydata
{


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

auto supported_layouts()
{
    using enum core::LayoutMode;
    std::vector layouts{AoSMapped};
    PHARE_WITH_MKN_GPU({
        layouts.emplace_back(AoSTS);
        layouts.emplace_back(AoSPCTS);
    })
    return layouts;
}



PYBIND11_MODULE(cpp_etc, m, py::mod_gil_not_used())
{
    auto samrai_restart_file = [](std::string path) {
        return PHARE::amr::HierarchyRestarter::getRestartFileFullPath(path);
    };

    py::class_<core::Span<double>, py::smart_holder>(m, "Span");
    m.def("makeSpan", makePySpan<double>);
    py::class_<PyArrayWrapper<double>, py::smart_holder, core::Span<double>>(m, "PyWrapper");


    m.def("mpi_size", []() { return mpi::size(); });
    m.def("mpi_rank", []() { return mpi::rank(); });
    m.def("mpi_barrier", []() { mpi::barrier(); });
    m.def("mpi_initialized", []() { return mpi::is_init(); });

    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, py::smart_holder>(m, "AMRHierarchy");
    m.def("make_hierarchy", []() { return PHARE::make_hierarchy(); });

    m.def("makePyArrayWrapper", makePyArrayWrapper<double>);

    m.def("phare_deps", []() {
        std::unordered_map<std::string, std::string> versions{{"pybind", pybind_version()},
                                                              {"samrai", samrai_version()}};
        _PHARE_WITH_HIGHFIVE(versions["highfive"] = HIGHFIVE_VERSION_STRING);
        return versions;
    });

    m.def("samrai_restart_file", samrai_restart_file);

    m.def("restart_path_for_time", [](std::string path, double timestamp) {
        return PHARE::amr::Hierarchy::restartFilePathForTime(path, timestamp);
    });

    m.def("phare_build_config", []() { return PHARE::build_config(); });

    m.def("patch_data_ids", [&](std::string const& path) -> std::vector<int> {
        _PHARE_WITH_HIGHFIVE({
            auto const& restart_file = samrai_restart_file(path);
            PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadOnly,
                                                 /*para=*/false};
            return h5File.read_data_set<int>("/phare/patch/ids");
        });

        throw std::runtime_error("PHARE not built with highfive support");
    });
    m.def("serialized_simulation_string", [&](std::string const& path) -> std::string {
        _PHARE_WITH_HIGHFIVE({
            auto const& restart_file = samrai_restart_file(path);
            PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadOnly,
                                                 /*para=*/false};
            return h5File.read_attribute("/phare", "serialized_simulation");
        });

        throw std::runtime_error("PHARE not built with highfive support");
    });

    using enum core::LayoutMode;
    py::native_enum<core::LayoutMode>(m, "LayoutMode", "enum.Enum")
        .value("AoSMapped", AoSMapped)
        .value("AoSTS", AoSTS)
        .value("AoSPCTS", AoSPCTS)
        .finalize();


    m.def("supported_layouts", supported_layouts);

    py::native_enum<AllocatorMode>(m, "AllocatorMode", "enum.Enum")
        .value("CPU", AllocatorMode::CPU)
        .value("GPU_UNIFIED", AllocatorMode::GPU_UNIFIED)
        .finalize();

    py::native_enum<MHDOpts::TimeIntegratorType>(m, "TimeIntegratorType", "enum.Enum")
        .value("Default", MHDOpts::TimeIntegratorType::Default)
        .value("Euler", MHDOpts::TimeIntegratorType::Euler)
        .value("TVDRK2", MHDOpts::TimeIntegratorType::TVDRK2)
        .value("TVDRK3", MHDOpts::TimeIntegratorType::TVDRK3)
        .value("SSPRK4_5", MHDOpts::TimeIntegratorType::SSPRK4_5)
        .finalize();

    py::native_enum<MHDOpts::ReconstructionType>(m, "ReconstructionType", "enum.Enum")
        .value("Default", MHDOpts::ReconstructionType::Default)
        .value("Constant", MHDOpts::ReconstructionType::Constant)
        .value("Linear", MHDOpts::ReconstructionType::Linear)
        .value("WENO3", MHDOpts::ReconstructionType::WENO3)
        .value("WENOZ", MHDOpts::ReconstructionType::WENOZ)
        .value("MP5", MHDOpts::ReconstructionType::MP5)
        .finalize();

    py::native_enum<MHDOpts::SlopeLimiterType>(m, "SlopeLimiterType", "enum.Enum")
        .value("None", MHDOpts::SlopeLimiterType::None) // accessed via getattr(), not dot syntax
        .value("VanLeer", MHDOpts::SlopeLimiterType::VanLeer)
        .value("MinMod", MHDOpts::SlopeLimiterType::MinMod)
        .finalize();

    py::native_enum<MHDOpts::RiemannSolverType>(m, "RiemannSolverType", "enum.Enum")
        .value("Default", MHDOpts::RiemannSolverType::Default)
        .value("Rusanov", MHDOpts::RiemannSolverType::Rusanov)
        .value("HLL", MHDOpts::RiemannSolverType::HLL)
        .value("HLLD", MHDOpts::RiemannSolverType::HLLD)
        .finalize();

    py::class_<SimOpts>(m, "SimOpts")
        .def(py::init<>())
        .def_readwrite("dimension", &SimOpts::dimension)
        .def_readwrite("interp_order", &SimOpts::interp_order)
        .def_readwrite("nbRefinedPart", &SimOpts::nbRefinedPart)
        .def_readwrite("layout_mode", &SimOpts::layout_mode)
        .def_readwrite("alloc_mode", &SimOpts::alloc_mode)
        .def_readwrite("time_integrator_type", &SimOpts::time_integrator_type)
        .def_readwrite("reconstruction_type", &SimOpts::reconstruction_type)
        .def_readwrite("slope_limiter_type", &SimOpts::slope_limiter_type)
        .def_readwrite("riemann_solver_type", &SimOpts::riemann_solver_type)
        .def_readwrite("Hall", &SimOpts::Hall)
        .def_readwrite("Resistivity", &SimOpts::Resistivity)
        .def_readwrite("HyperResistivity", &SimOpts::HyperResistivity)
        .def("__eq__", [](SimOpts const& a, SimOpts const& b) { return a == b; });
}

} // namespace PHARE::pydata
