#include "samrai.hpp"


namespace PHARE
{

SamraiLifeCycle::SamraiLifeCycle(int argc, char** argv)
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
            phlop::ScopeTimerMan::INSTANCE()
                .file_name(".phare/timings/rank." + std::to_string(core::mpi::rank()) + ".txt")
                .init(); //
    )
}

SamraiLifeCycle::~SamraiLifeCycle()
{
    PHARE_WITH_PHLOP(phlop::ScopeTimerMan::reset());
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
}

void SamraiLifeCycle::reset()
{
    PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::startup();
    SAMRAI::tbox::RestartManager::getManager()->clearRestartItems();
}




} // namespace PHARE
