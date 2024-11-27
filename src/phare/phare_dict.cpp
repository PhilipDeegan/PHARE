//
//  This exe is provided to be able to run PHARE from restart files without python
//

#include "phare/phare.hpp"
#include "simulator/simulator.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "initializer/data_provider.hpp"

#include <atomic>
#include <csignal>

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
        throw std::runtime_error("PHARE Error: Failed to register SIGINT signal handler");
    if (std::signal(SIGABRT, signal_handler) == SIG_ERR)
        throw std::runtime_error("PHARE Error: Failed to register SIGABRT signal handler");
    if (argc != 2)
        throw std::runtime_error("PHARE Error: Executable expects one argument");
    std::string const welcome = R"~(
                  _____   _    _            _____   ______
                 |  __ \ | |  | |    /\    |  __ \ |  ____|
                 | |__) || |__| |   /  \   | |__) || |__
                 |  ___/ |  __  |  / /\ \  |  _  / |  __|
                 | |     | |  | | / ____ \ | | \ \ | |____
                 |_|     |_|  |_|/_/    \_\|_|  \_\|______|)~";
    std::cout << welcome << "\n\n";
    PHARE::SamraiLifeCycle slc{argc, argv};
    PHARE::initializer::PHAREDictHandler::load(argv[1]);
    auto hierarchy = PHARE::amr::Hierarchy::make();
    auto simulator = PHARE::getSimulator(hierarchy);
    std::cout << PHARE::core::to_str(*simulator) << "\n";
    simulator->initialize();
    PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
    while (simulator->currentTime() < simulator->endTime())
    {
        if (gSignalStatus)
            return gSignalStatus;

        simulator->dump(simulator->currentTime(), simulator->timeStep());
        simulator->advance(simulator->timeStep());
    }
    return gSignalStatus;
}
