
#ifndef PHARE_PHARE_INCLUDE_HPP
#define PHARE_PHARE_INCLUDE_HPP

#include <memory>
#include <cassert>
#include <iostream>

#include "core/def/phlop.hpp" // scope timing
#include "core/utilities/algorithm.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "simulator/simulator.hpp"
#include "initializer/python_data_provider.hpp"


namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(std::string const& message, std::string const& filename, const int line)
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
    SamraiLifeCycle(int argc = 0, char** argv = nullptr)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();
        // uncomment next line for debugging samrai issues
        // SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort();
        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
        PHARE_WITH_PHLOP( //
            if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
                phlop::threaded::ScopeTimerMan::INSTANCE()
                    .file_name(".phare/timings/rank." + std::to_string(core::mpi::rank()) + ".txt")
                    .init(); //
        )
    }

    ~SamraiLifeCycle()
    {
        PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
    }

    static void reset()
    {
        PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
};



std::unique_ptr<PHARE::initializer::DataProvider> fromCommandLine(int argc, char** argv)
{
    using dataProvider [[maybe_unused]] = std::unique_ptr<PHARE::initializer::DataProvider>;

    switch (argc)
    {
        case 1: return nullptr;
        case 2:
            std::string arg = argv[1];
            auto moduleName = arg.substr(0, arg.find_last_of("."));
            if (arg.substr(arg.find_last_of(".") + 1) == "py")
            {
                std::replace(moduleName.begin(), moduleName.end(), '/', '.');
                std::cout << "python input detected, building with python provider...\n";
                return std::make_unique<PHARE::initializer::PythonDataProvider>(moduleName);
            }

            break;
    }
    return nullptr;
}

} // namespace PHARE


#endif /*PHARE_PHARE_INCLUDE_H*/
