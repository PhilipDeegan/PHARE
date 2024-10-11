// Full simulator test - with custom init - and no advances!


struct PHARE_Testables;
#define PHARE_FRIEND_CLASS_HACKERY friend struct ::PHARE_Testables
#include "phare/phare.hpp"
#undef PHARE_FRIEND_CLASS_HACKERY

#include "gtest/gtest.h"
#include "tests/initializer/init_functions.hpp"
#include "initializer/python_data_provider.hpp"


using namespace PHARE::core;

int nbrPartPerCell = 1000;

namespace PHARE::initializer
{

template<std::size_t dim = 2>
PHARE::initializer::PHAREDict createDict()
{
    PHAREDictHandler::INSTANCE().init();

    char const* name = "job";
    PythonDataProvider pydp{name};
    pydp.read();
    return PHAREDictHandler::INSTANCE().dict();
}

} // namespace PHARE::initializer



struct PHARE_Testables
{
    std::size_t constexpr static dim          = 2;
    std::size_t constexpr static interp       = 1;
    std::size_t constexpr static nbRefinePart = 4;
    using MultiPhysicsIntegrator_t =
        typename PHARE::solver::PHARE_Types<dim, interp, nbRefinePart>::MultiPhysicsIntegrator;


    void operator()()
    {
        sim.initialize();
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();

        auto& integrator = sim.integrator_;
        auto& timeRefInt = integrator->timeRefIntegrator_;

        auto& multi = *dynamic_cast<MultiPhysicsIntegrator_t*>(timeRefInt.get());
        PHARE_LOG_LINE_SS(multi.nbrOfLevels_);
    }

    PHARE::initializer::PHAREDict dict               = PHARE::initializer::createDict();
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy = PHARE::amr::Hierarchy::make();
    PHARE::Simulator<dim, interp, nbRefinePart> sim{dict, hierarchy};
};

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle slc{argc, argv};
    PHARE_Testables{}();
    return 0;
}
