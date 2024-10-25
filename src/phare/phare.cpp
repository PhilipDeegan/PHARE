#include "phare/phare.hpp"
#include "simulator/simulator.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "initializer/python_data_provider.hpp"

#include "core/logger.hpp"

#include <algorithm>
#include <csignal>
#include <atomic>

namespace PHARE
{
std::unique_ptr<PHARE::ISimulator> getSimulator(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
{
    PHARE::initializer::PHAREDict const& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim           = theDict["simulation"]["dimension"].template to<int>();
    auto interpOrder   = theDict["simulation"]["interp_order"].template to<int>();
    auto nbRefinedPart = theDict["simulation"]["refined_particle_nbr"].template to<int>();

    return core::makeAtRuntime<SimulatorMaker>(dim, interpOrder, nbRefinedPart,
                                               SimulatorMaker{hierarchy});
}
} /* namespace PHARE */


namespace
{
std::atomic<int> gSignalStatus = 0;
}

void signal_handler(int signal)
{
    gSignalStatus = signal;
}



int main(int argc, char** argv)
{
    if (std::signal(SIGINT, signal_handler) == SIG_ERR)
    {
        throw std::runtime_error("PHARE Error: Failed to register SIGINT signal handler");
    }
    if (std::signal(SIGABRT, signal_handler) == SIG_ERR)
    {
        throw std::runtime_error("PHARE Error: Failed to register SIGABRT signal handler");
    }

    std::string const welcome = R"~(
                  _____   _    _            _____   ______
                 |  __ \ | |  | |    /\    |  __ \ |  ____|
                 | |__) || |__| |   /  \   | |__) || |__
                 |  ___/ |  __  |  / /\ \  |  _  / |  __|
                 | |     | |  | | / ____ \ | | \ \ | |____
                 |_|     |_|  |_|/_/    \_\|_|  \_\|______|)~";
    std::cout << welcome;
    std::cout << "\n";
    std::cout << "\n";

    PHARE::SamraiLifeCycle slc{argc, argv};

    std::cerr << "creating python data provider\n";
    auto provider = PHARE::fromCommandLine(argc, argv);

    std::cerr << "reading user inputs...";
    provider->read();
    std::cerr << "done!\n";

    auto& dictHandler = PHARE::initializer::PHAREDictHandler::INSTANCE();

    auto hierarchy = PHARE::amr::Hierarchy::make();

    auto simulator = PHARE::getSimulator(hierarchy);

    std::cout << PHARE::core::to_str(*simulator) << "\n";

    simulator->initialize();

    dictHandler.stop();
    provider.release();

    [[maybe_unused]] auto time = simulator->startTime();

    while (simulator->currentTime() < simulator->endTime())
    {
        if (gSignalStatus)
            return gSignalStatus;

        simulator->dump(simulator->currentTime(), simulator->timeStep());
        simulator->advance(simulator->timeStep());
    }

    return gSignalStatus;
}