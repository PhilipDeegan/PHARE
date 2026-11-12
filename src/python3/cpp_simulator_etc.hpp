#ifndef PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP

#include "python3/cpp_simulator.hpp"

#include "python3/particles.hpp"     // IWYU pragma: keep
#include "python3/patch_level.hpp"   // IWYU pragma: keep
#include "python3/data_wrangler.hpp" // IWYU pragma: keep

namespace PHARE::pydata
{

template<auto opts>
void declare_etc(py::module& m)
{
    using Sim = Simulator<opts>;

    using DW         = DataWrangler<opts>;
    std::string name = "DataWrangler";

    auto dwCls
        = py::class_<DW, py::smart_holder>(m, name.c_str(), py::module_local())
              .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
              .def("sync", &DW::sync)
              .def("getNumberOfLevels", &DW::getNumberOfLevels);

    if constexpr (has_hybrid_v<opts>)
        dwCls.def("getHybridPatchLevel", &DW::getHybridPatchLevel);
    if constexpr (has_mhd_v<opts>)
        dwCls.def("getMHDPatchLevel", &DW::getMHDPatchLevel);


    if constexpr (has_hybrid_v<opts>)
    {
        using HybPL = PatchLevel<typename solver::PHARE_Types<opts>::Hybrid::Model_t>;
        name        = "HybridPatchLevel";
        if constexpr (defaultNbrRefinedParts(opts.dimension, opts.interp_order)
                      == opts.nbRefinedPart) // register once!
            py::class_<HybPL, py::smart_holder>(m, name.c_str(), py::module_local())
                .def("getB", &HybPL::getB, py::arg("component"))
                .def("getE", &HybPL::getE, py::arg("component"))
                .def("getVi", &HybPL::getVi, py::arg("component"))
                .def("getN", &HybPL::getN, py::arg("pop_name"))
                .def("getNi", &HybPL::getNi)
                .def("getFlux", &HybPL::getFlux, py::arg("component"), py::arg("pop_name"))
                .def("getParticles", &HybPL::getParticles, py::arg("pop_name"));

        using _Splitter = PHARE::amr::Splitter<core::DimConst<opts.dimension>,
                                               core::InterpConst<opts.interp_order>,
                                               core::RefinedParticlesConst<opts.nbRefinedPart>>;
        name            = "Splitter";

        py::class_<_Splitter, py::smart_holder>(m, name.c_str(), py::module_local())
            .def(py::init<>())
            .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
            .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

        name = "split_pyarray_particles";
        m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
    }

    if constexpr (has_mhd_v<opts>)
    {
        using MHDPL = PatchLevel<typename solver::PHARE_Types<opts>::MHD::Model_t>;
        name        = "MHDPatchLevel";
        py::class_<MHDPL, py::smart_holder>(m, name.c_str(), py::module_local())
            .def("getRho", &MHDPL::getRho)
            .def("getP", &MHDPL::getP)
            .def("getEtot", &MHDPL::getEtot)
            .def("getV", &MHDPL::getV, py::arg("component"))
            .def("getB", &MHDPL::getB, py::arg("component"))
            .def("getRhoV", &MHDPL::getRhoV, py::arg("component"))
            .def("getE", &MHDPL::getE, py::arg("component"))
            .def("getJ", &MHDPL::getJ, py::arg("component"));
    }
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP*/
