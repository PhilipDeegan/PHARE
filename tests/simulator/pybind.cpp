
#include "phare/include.h"

#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(const std::string& message, const std::string& filename, const int line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle
{
public:
    static SamraiLifeCycle& INSTANCE()
    {
        static SamraiLifeCycle i;
        return i;
    }

    SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAI_MPI::init(0, nullptr);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();

        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    }
    ~SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    // the simulator must be destructed before the variable database is reset, or segfault.
    void reset()
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
};


template<size_t dim, size_t interp>
void declareParticlesDim(py::module& m)
{
    using T          = ContiguousParticles<dim, interp>;
    std::string name = "ContiguousParticles_" + std::to_string(dim) + "_" + std::to_string(interp);
    py::class_<T, std::shared_ptr<T>>(m, name.c_str())
        .def(py::init<size_t>())
        .def_readwrite("iCell", &T::iCell)
        .def_readwrite("delta", &T::delta)
        .def_readwrite("weight", &T::weight)
        .def_readwrite("charge", &T::charge)
        .def_readwrite("v", &T::v)
        .def("size", &T::size)
        .def("split", &T::split);
}

template<size_t interp>
void declareParticles(py::module& m)
{
    declareParticlesDim<1, interp>(m);
    // declareParticlesDim<2, interp>(m);
    // declareParticlesDim<3, interp>(m);
}

PYBIND11_MODULE(test_simulator, m)
{
    SamraiLifeCycle::INSTANCE(); // init
    declareParticles<1>(m);
    declareParticles<2>(m);
    declareParticles<3>(m);

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");

    py::class_<ISimulator, std::shared_ptr<ISimulator>>(m, "ISimulator")
        .def("initialize", &PHARE::ISimulator::initialize)
        .def("advance", &PHARE::ISimulator::advance)
        .def("startTime", &PHARE::ISimulator::startTime)
        .def("currentTime", &PHARE::ISimulator::currentTime)
        .def("endTime", &PHARE::ISimulator::endTime)
        .def("timeStep", &PHARE::ISimulator::timeStep)
        .def("to_str", &PHARE::ISimulator::to_str);

    py::class_<RuntimeDiagnosticInterface, std::shared_ptr<RuntimeDiagnosticInterface>>(
        m, "IDiagnosticsManager")
        .def("dump", &RuntimeDiagnosticInterface::dump);

    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });
    m.def("make_simulator", [](std::shared_ptr<PHARE::amr::Hierarchy>& hier) {
        auto sim = PHARE::getSimulator(hier);
        auto ptr = sim.get();
        sim.release();
        return std::shared_ptr<ISimulator>{ptr};
    });
    m.def("make_diagnostic_manager", [](std::shared_ptr<ISimulator> const& sim,
                                        std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::make_shared<RuntimeDiagnosticInterface>(*sim, *hier);
    });

    m.def("mpi_size", []() {
        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        return mpi_size;
    });
    m.def("mpi_rank", []() {
        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        return mpi_rank;
    });
    m.def("reset", []() {
        py::gil_scoped_release release;
        SamraiLifeCycle::INSTANCE().reset();
    });
}
} // namespace PHARE
