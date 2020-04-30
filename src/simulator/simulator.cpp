#include "simulator.h"


namespace PHARE
{
/*
    template definitions for library for methods in this file
    useful for minimizing compile times during dev
*/

template<std::size_t dim, std::size_t io, size_t nb>
Simulator<dim, io, nb>::Simulator(PHARE::initializer::PHAREDict dict,
                                  std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
    : hierarchy_{hierarchy}
    , modelNames_{"HybridModel"}
    , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
    , messengerFactory_{descriptors_}
    , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
    , dt_{dict["simulation"]["time_step"].template to<double>()}
    , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
    , finalTime_{dt_ * timeStepNbr_}
    , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(
          dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>())}
{
    if (find_model("HybridModel"))
    {
        hybridModel_ = std::make_shared<HybridModel>(
            dict["simulation"], std::make_shared<typename HybridModel::resources_manager_type>());


        hybridModel_->resourcesManager->registerResources(hybridModel_->state);

        // we register the hybrid model for all possible levels in the hierarchy
        // since for now it is the only model available
        // same for the solver
        multiphysInteg_->registerModel(0, maxLevelNumber_ - 1, hybridModel_);
        multiphysInteg_->registerAndInitSolver(
            0, maxLevelNumber_ - 1, std::make_unique<SolverPPC>(dict["simulation"]["solverPPC"]));
        multiphysInteg_->registerAndSetupMessengers(messengerFactory_);


        auto startTime = 0.; // TODO make it runtime
        auto endTime   = 0.; // TODO make it runtime


        integrator_ = std::make_unique<PHARE::amr::DimIntegrator<dimension>>(
            dict, hierarchy, multiphysInteg_, multiphysInteg_, startTime, endTime);
    }
    else
        throw std::runtime_error("unsupported model");
}



template<std::size_t dim, std::size_t io, size_t nb>
std::string Simulator<dim, io, nb>::to_str()
{
    std::stringstream ss;
    ss << "PHARE SIMULATOR\n";
    ss << "------------------------------------\n";
    ss << "interpolation order  : " << interp_order << "\n";
    ss << "dimension            : " << dimension << "\n";
    ss << "time step            : " << dt_ << "\n";
    ss << "number of time steps : " << timeStepNbr_ << "\n";
    ss << core::to_str(hybridModel_->state);
    return ss.str();
}



template<std::size_t dim, std::size_t io, size_t nb>
void Simulator<dim, io, nb>::initialize()
{
    if (integrator_ != nullptr)
        integrator_->initialize();
    else
        throw std::runtime_error("Error - Simulator has no integrator");
}



template<std::size_t dim, std::size_t io, size_t nb>
double Simulator<dim, io, nb>::startTime()
{
    return 0.;
}
template<std::size_t dim, std::size_t io, size_t nb>
double Simulator<dim, io, nb>::endTime()
{
    return finalTime_;
}
template<std::size_t dim, std::size_t io, size_t nb>
double Simulator<dim, io, nb>::timeStep()
{
    return dt_;
}
template<std::size_t dim, std::size_t io, size_t nb>
double Simulator<dim, io, nb>::currentTime()
{
    return currentTime_;
}


template<std::size_t dim, std::size_t io, size_t nb>
void Simulator<dim, io, nb>::advance()
{
    currentTime_ += dt_;
    timeStepNbr_++;
}


template class Simulator<1, 1, 2>;
template class Simulator<1, 2, 2>;
template class Simulator<1, 3, 2>;


} /* namespace PHARE */
