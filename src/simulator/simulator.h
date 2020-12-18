
#ifndef PHARE_SIMULATOR_SIMULATOR_H
#define PHARE_SIMULATOR_SIMULATOR_H

#include "include.h"

#include "phare_core.h"
#include "phare_types.h"

#include "amr/tagging/tagger_factory.h"
#include <chrono>


namespace PHARE
{
template<typename Float>
class ISimulator
{
public:
    virtual Float startTime()   = 0;
    virtual Float endTime()     = 0;
    virtual Float currentTime() = 0;
    virtual Float timeStep()    = 0;

    virtual void initialize()       = 0;
    virtual Float advance(Float dt) = 0;

    virtual std::vector<int> const& domainBox() const   = 0;
    virtual std::vector<Float> const& cellWidth() const = 0;
    virtual std::size_t interporder() const             = 0;

    virtual std::string to_str() = 0;

    virtual ~ISimulator() {}
    virtual void dump(double timestamp, double timestep) {} // overriding optional
};

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
class Simulator : public ISimulator<Float_>
{
public:
    using Float = Float_;

    Float startTime() override { return 0.; }
    Float endTime() override { return finalTime_; }
    Float timeStep() override { return dt_; }
    Float currentTime() override { return currentTime_; }

    void initialize() override;
    Float advance(Float dt) override;

    std::vector<int> const& domainBox() const override { return hierarchy_->domainBox(); }
    std::vector<Float> const& cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    auto& getHybridModel() { return hybridModel_; }
    auto& getMHDModel() { return mhdModel_; }
    auto& getMultiPhysicsIntegrator() { return multiphysInteg_; }

    std::string to_str() override;

    void dump(double timestamp, double timestep) override { dMan->dump(timestamp, timestep); }

    Simulator(PHARE::initializer::PHAREDict dict,
              std::shared_ptr<PHARE::amr::Hierarchy<Float>> const& hierarchy);

    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using SAMRAITypes = PHARE::amr::SAMRAI_Types;
    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart, Float>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<SAMRAITypes, Float>;
    using HybridModel    = typename PHARETypes::HybridModel_t;
    using MHDModel       = typename PHARETypes::MHDModel_t;

    using SolverMHD = typename PHARETypes::SolverMHD_t;
    using SolverPPC = typename PHARETypes::SolverPPC_t;

    using MessengerFactory       = typename PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;

    using SimFunctorParams = typename core::PHARE_Sim_Types::SimFunctorParams;
    using SimFunctors      = typename core::PHARE_Sim_Types::SimulationFunctors;

    using Integrator  = PHARE::amr::Integrator<Float, dimension>;
    using Hierarchy_t = PHARE::amr::Hierarchy<Float>;

private:
    auto find_model(std::string name);

    std::shared_ptr<Hierarchy_t> hierarchy_;
    std::unique_ptr<Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];
    int maxLevelNumber_;
    Float dt_;
    int timeStepNbr_           = 0;
    Float finalTime_           = 0;
    Float currentTime_         = 0;
    bool isInitialized         = false;
    std::size_t fineDumpLvlMax = 0;

    // physical models that can be used
    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;

    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;

    SimFunctors functors_;

    SimFunctors functors_setup(PHARE::initializer::PHAREDict const& dict)
    {
        return {{"pre_advance", {/*empty vector*/}}};
    }

    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_{nullptr};
};


//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
Simulator<_dimension, _interp_order, _nbRefinedPart, Float_>::Simulator(
    PHARE::initializer::PHAREDict dict, std::shared_ptr<Hierarchy_t> const& hierarchy)
    : hierarchy_{hierarchy}
    , modelNames_{"HybridModel"}
    , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
    , messengerFactory_{descriptors_}
    , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
    , dt_{dict["simulation"]["time_step"].template to<Float_>()}
    , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
    , finalTime_{dt_ * timeStepNbr_}
    , functors_{functors_setup(dict)}
    , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(
          dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>(), functors_)}
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
            0, maxLevelNumber_ - 1, std::make_unique<SolverPPC>(dict["simulation"]["algo"]));

        multiphysInteg_->registerAndSetupMessengers(messengerFactory_);

        // hard coded for now, should get some params later from the dict
        auto hybridTagger_ = amr::TaggerFactory<PHARETypes>::make("HybridModel", "default");
        multiphysInteg_->registerTagger(0, maxLevelNumber_ - 1, std::move(hybridTagger_));


        auto startTime = 0.; // TODO make it runtime
        auto endTime   = 0.; // TODO make it runtime


        integrator_ = std::make_unique<Integrator>(dict, hierarchy, multiphysInteg_,
                                                   multiphysInteg_, startTime, endTime);

        if (dict["simulation"].contains("diagnostics"))

        {
            auto& diagDict = dict["simulation"]["diagnostics"];

            dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(
                *hierarchy_, *hybridModel_, diagDict);

            if (diagDict.contains("fine_dump_lvl_max"))
            {
                auto fine_dump_lvl_max = diagDict["fine_dump_lvl_max"].template to<int>();

                if (fine_dump_lvl_max > 0)
                { // copy for later
                    this->fineDumpLvlMax = static_cast<std::size_t>(fine_dump_lvl_max);
                    functors_["pre_advance"]["fine_dump"] = [&](SimFunctorParams const& params) {
                        std::size_t level_nbr = params["level_nbr"].template to<int>();
                        auto timestamp        = params["timestamp"].template to<double>();

                        if (this->fineDumpLvlMax >= level_nbr)
                            this->dMan->dump_level(level_nbr, timestamp);
                    };
                }
            }
        }
    }
    else
        throw std::runtime_error("unsupported model");
}



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
std::string Simulator<_dimension, _interp_order, _nbRefinedPart, Float_>::to_str()
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




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
void Simulator<_dimension, _interp_order, _nbRefinedPart, Float_>::initialize()
{
    try
    {
        if (isInitialized)
            std::runtime_error("cannot initialize  - simulator already isInitialized");

        if (integrator_ != nullptr)
            integrator_->initialize();
        else
            throw std::runtime_error("Error - Simulator has no integrator");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    isInitialized = true;
}



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
Float_ Simulator<_dimension, _interp_order, _nbRefinedPart, Float_>::advance(Float_ dt)
{
    try
    {
        if (integrator_)
        {
            auto dt_new = integrator_->advance(dt);
            currentTime_ += dt;
            return dt_new;
        }
        else
            throw std::runtime_error("Error - no valid integrator in the simulator");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }
}



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         typename Float_>
auto Simulator<_dimension, _interp_order, _nbRefinedPart, Float_>::find_model(std::string name)
{
    return std::find(std::begin(modelNames_), std::end(modelNames_), name) != std::end(modelNames_);
}



template<typename Float>
struct SimulatorMaker
{
    using Hierarchy_t = PHARE::amr::Hierarchy<Float>;

    SimulatorMaker(std::shared_ptr<Hierarchy_t>& hierarchy)
        : hierarchy_{hierarchy}
    {
    }

    std::shared_ptr<Hierarchy_t>& hierarchy_;

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    std::unique_ptr<ISimulator<Float>>
    operator()(std::size_t userDim, std::size_t userInterpOrder, std::size_t userNbRefinedPart,
               Dimension dimension, InterpOrder interp_order, NbRefinedPart nbRefinedPart)
    {
        if (userDim == dimension() and userInterpOrder == interp_order()
            and userNbRefinedPart == nbRefinedPart())
        {
            std::size_t constexpr d  = dimension();
            std::size_t constexpr io = interp_order();
            std::size_t constexpr nb = nbRefinedPart();

            PHARE::initializer::PHAREDict& theDict
                = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            return std::make_unique<Simulator<d, io, nb, Float>>(theDict, hierarchy_);
        }
        else
        {
            return nullptr;
        }
    }
};


template<std::size_t dimension, std::size_t interp_order, std::size_t nbRefinedPart, typename Float>
class SimulatorCaster
{
public:
    using Simulator_t = Simulator<dimension, interp_order, nbRefinedPart, Float>;

    SimulatorCaster(std::shared_ptr<ISimulator<Float>> const& _simulator)
        : simulator{_simulator}
    {
    }

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    Simulator_t* operator()(std::size_t userDim, std::size_t userInterpOrder,
                            std::size_t userNbRefinedPart, Dimension dimension_fn,
                            InterpOrder interp_order_fn, NbRefinedPart nbRefinedPart_fn)
    {
        if (userDim == dimension_fn() and userInterpOrder == interp_order_fn()
            and userNbRefinedPart == nbRefinedPart_fn())
        {
            std::size_t constexpr d  = dimension_fn();
            std::size_t constexpr io = interp_order_fn();
            std::size_t constexpr nb = nbRefinedPart_fn();

            // extra if constexpr as cast is templated and not generic interface
            if constexpr (d == dimension and io == interp_order and nb == nbRefinedPart)
                return dynamic_cast<Simulator_t*>(simulator.get());
        }
        return nullptr;
    }

private:
    std::shared_ptr<ISimulator<Float>> const& simulator;
};

template<typename Float>
std::unique_ptr<PHARE::ISimulator<Float>>
getSimulator(std::shared_ptr<PHARE::amr::Hierarchy<Float>>& hierarchy)
{
    PHARE::initializer::PHAREDict& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim           = theDict["simulation"]["dimension"].template to<int>();
    auto interpOrder   = theDict["simulation"]["interp_order"].template to<int>();
    auto nbRefinedPart = theDict["simulation"]["refined_particle_nbr"].template to<int>();

    return core::makeAtRuntime<Float, SimulatorMaker<Float>>(dim, interpOrder, nbRefinedPart,
                                                             SimulatorMaker<Float>{hierarchy});
}


template<std::size_t dim, std::size_t interp, std::size_t nbRefinedPart, typename Float>
std::unique_ptr<Simulator<dim, interp, nbRefinedPart, Float>>
makeSimulator(std::shared_ptr<amr::Hierarchy<Float>> const& hierarchy)
{
    return std::make_unique<Simulator<dim, interp, nbRefinedPart, Float>>(
        initializer::PHAREDictHandler::INSTANCE().dict(), hierarchy);
}


} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_H*/
