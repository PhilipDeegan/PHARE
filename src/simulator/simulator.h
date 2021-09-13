
#ifndef PHARE_SIMULATOR_SIMULATOR_H
#define PHARE_SIMULATOR_SIMULATOR_H

#include "include.h"

#include "phare_core.h"
#include "phare_types.h"

#include "core/utilities/types.h"
#include "core/utilities/mpi_utils.h"
#include "core/utilities/timestamps.h"
#include "amr/tagging/tagger_factory.h"

#include <chrono>
#include <exception>


#if defined(HAVE_UMPIRE)
#include "umpire/ResourceManager.hpp"
#endif


namespace PHARE{
auto static_init(){
#if defined(HAVE_UMPIRE)

	
#endif  
  return true;
}
}


namespace PHARE
{
class ISimulator
{
public:
    virtual double startTime()   = 0;
    virtual double endTime()     = 0;
    virtual double currentTime() = 0;
    virtual double timeStep()    = 0;

    virtual void initialize()         = 0;
    virtual double advance(double dt) = 0;

    virtual std::vector<int> const& domainBox() const    = 0;
    virtual std::vector<double> const& cellWidth() const = 0;
    virtual std::size_t interporder() const              = 0;

    virtual std::string to_str() = 0;

    virtual ~ISimulator() {}
    virtual bool dump(double timestamp, double timestep) { return false; } // overriding optional
};

template<typename StateType>
struct State
{
    static constexpr std::size_t dimension = StateType::dimension;
    using HybridModel                      = typename StateType::HybridModel;
    using MultiPhysicsIntegrator           = typename StateType::MultiPhysicsIntegrator;
    using SolverPPC                        = typename StateType::SolverPPC;
    using Integrator                       = PHARE::amr::Integrator<dimension>;

    template<typename SimFunctors>
    State(PHARE::initializer::PHAREDict const& dict, SimFunctors& functors)
        : multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(dict["simulation"], functors)}
    {
    }

    template<typename OtherType>
    State(State<OtherType>& that)
        : hybridModel_{that.hybridModel_}
        , multiphysInteg_{that.multiphysInteg_}
        , integrator_{that.integrator_}
    {
    }

    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_;
    std::shared_ptr<Integrator> integrator_;
};


// template<std::size_t _dimension, std::size_t interp_order, std::size_t nbRefinedPart>
// struct InitState // always CPU
// {
//     static constexpr std::size_t dimension = _dimension;
//     using PHARETypes             = solver::CPU_Types<dimension, interp_order, nbRefinedPart>;
//     using HybridModel            = typename PHARETypes::HybridModel_t;
//     using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;
//     using SolverPPC              = typename PHARETypes::SolverPPC_t;
// };

template<std::size_t _dimension, std::size_t interp_order, std::size_t nbRefinedPart, bool offload>
struct WorkingState // might be CPU
{
    static constexpr std::size_t dimension = _dimension;
    using PHARETypes             = PHARE_Types<dimension, interp_order, nbRefinedPart, offload>;
    using HybridModel            = typename PHARETypes::HybridModel_t;
    using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;
    using SolverPPC              = typename PHARETypes::SolverPPC_t;

    // static constexpr bool is_cpu
    //     = std::is_same_v<HybridModel,
    //                      typename InitState<dimension, interp_order,
    //                      nbRefinedPart>::HybridModel>;
};

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         bool offload = false>
class Simulator : public ISimulator
{
public:
    double startTime() override { return 0.; }
    double endTime() override { return finalTime_; }
    double timeStep() override { return dt_; }
    double currentTime() override { return currentTime_; }

    void initialize() override;
    double advance(double dt) override;

    std::vector<int> const& domainBox() const override { return hierarchy_->domainBox(); }
    std::vector<double> const& cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    auto& getHybridModel() { return state_->hybridModel_; }
    auto& getMHDModel() { return mhdModel_; }
    auto& getMultiPhysicsIntegrator() { return state_->multiphysInteg_; }

    std::string to_str() override;

    bool dump(double timestamp, double timestep) override
    {
        // if constexpr (!WorkingState_t::is_cpu) // load cpu from gpu for diag dumps
        //     solver::fill_state_from(state_->hybridModel_, hostState_->hybridModel_);
        return dMan->dump(timestamp, timestep);
    }

    Simulator(PHARE::initializer::PHAREDict const& dict,
              std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy);
    ~Simulator()
    {
        if (coutbuf != nullptr)
            std::cout.rdbuf(coutbuf);
    }

    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using PHARETypes = PHARE_Types<dimension, interp_order, nbRefinedPart, offload>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel    = typename PHARETypes::HybridModel_t;
    using MHDModel       = typename PHARETypes::MHDModel_t;

    using SolverMHD = typename PHARETypes::SolverMHD_t;
    using SolverPPC = typename PHARETypes::SolverPPC_t;

    using MessengerFactory       = typename PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;

    using SimFunctorParams = typename core::PHARE_Sim_Types::SimFunctorParams;
    using SimFunctors      = typename core::PHARE_Sim_Types::SimulationFunctors;

    // using InitState_t    = InitState<dimension, interp_order, nbRefinedPart>;
    using WorkingState_t = WorkingState<dimension, interp_order, nbRefinedPart, offload>;


private:
    auto find_model(std::string name);

    std::ofstream log_out{".log/" + std::to_string(core::mpi::rank()) + ".out"};
    std::streambuf* coutbuf = nullptr;
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];
    int maxLevelNumber_;
    double dt_;
    int timeStepNbr_           = 0;
    double finalTime_          = 0;
    double currentTime_        = 0;
    bool isInitialized         = false;
    std::size_t fineDumpLvlMax = 0;

    // physical models that can be used
    // std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;

    std::unique_ptr<PHARE::core::ITimeStamper> timeStamper;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;

    SimFunctors functors_;

    SimFunctors functors_setup(PHARE::initializer::PHAREDict const& dict)
    {
        return {{"pre_advance", {/*empty vector*/}}};
    }


    // std::shared_ptr<State<InitState_t>> hostState_;
    std::shared_ptr<State<WorkingState_t>> state_;
};

namespace
{
    inline auto logging(std::ofstream& log_out)
    {
        std::streambuf* buf = nullptr;
        if (std::optional<std::string> log = core::get_env("PHARE_LOG");
            log and *log == "RANK_FILES")
        {
            buf = std::cout.rdbuf();
            std::cout.rdbuf(log_out.rdbuf());
        }
        return buf;
    }
} // namespace
//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------

template<typename State_t, typename Hierarchy_t, typename MessengerFactory_t>
auto init_state(std::shared_ptr<State_t>&& state, Hierarchy_t& hierarchy,
                MessengerFactory_t& messengerFactory, PHARE::initializer::PHAREDict const& dict,
                int maxLevelNumber, std::string prefix = "")
{
    using HybridModel = typename State_t::HybridModel;
    using SolverPPC   = typename State_t::SolverPPC;

    auto& hybridModel_    = state->hybridModel_;
    auto& multiphysInteg_ = state->multiphysInteg_;
    auto& integrator_     = state->integrator_;

    hybridModel_ = std::make_shared<HybridModel>(
        dict["simulation"], std::make_shared<typename HybridModel::resources_manager_type>(prefix));
    hybridModel_->resourcesManager->registerResources(hybridModel_->state);

    // we register the hybrid model for all possible levels in the hierarchy
    // since for now it is the only model available same for the solver
    multiphysInteg_->registerModel(0, maxLevelNumber - 1, hybridModel_);
    multiphysInteg_->registerAndInitSolver(0, maxLevelNumber - 1,
                                           std::make_unique<SolverPPC>(dict["simulation"]["algo"]));
    multiphysInteg_->registerAndSetupMessengers(messengerFactory);

    auto startTime = 0.; // TODO make it runtime
    auto endTime   = 0.; // TODO make it runtime

    integrator_ = std::make_shared<typename State_t::Integrator>(
        dict, hierarchy, multiphysInteg_, multiphysInteg_, startTime, endTime);

    return state;
}

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         bool offload>
Simulator<_dimension, _interp_order, _nbRefinedPart, offload>::Simulator(
    PHARE::initializer::PHAREDict const& dict,
    std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
    : coutbuf{logging(log_out)}
    , hierarchy_{hierarchy}
    , modelNames_{"HybridModel"}
    , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
    , messengerFactory_{descriptors_}
    , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
    , dt_{dict["simulation"]["time_step"].template to<double>()}
    , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
    , finalTime_{dt_ * timeStepNbr_}
    , functors_{functors_setup(dict)}
{
    if (find_model("HybridModel"))
    {
        // hostState_ = init_state(std::make_shared<State<InitState_t>>(dict, functors_),
        // hierarchy_,
        //                         messengerFactory_, dict, maxLevelNumber_);

        // if constexpr (!WorkingState_t::is_cpu)
        state_ = init_state(std::make_shared<State<WorkingState_t>>(dict, functors_), hierarchy_,
                            messengerFactory_, dict, maxLevelNumber_);
        // else
        // state_ = std::make_shared<State<WorkingState_t>>(*hostState_);

        // hard coded for now, should get some params later from the dict
        auto hybridTagger_ = amr::TaggerFactory<PHARETypes>::make("HybridModel", "default");
        state_->multiphysInteg_->registerTagger(0, maxLevelNumber_ - 1, std::move(hybridTagger_));

        timeStamper = core::TimeStamperFactory::create(dict["simulation"]);

        if (dict["simulation"].contains("diagnostics"))
        {
            auto& diagDict = dict["simulation"]["diagnostics"];

            dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(
                *hierarchy_, *state_->hybridModel_, diagDict);

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
         bool offload>
std::string Simulator<_dimension, _interp_order, _nbRefinedPart, offload>::to_str()
{
    std::stringstream ss;
    ss << "PHARE SIMULATOR\n";
    ss << "------------------------------------\n";
    ss << "interpolation order  : " << interp_order << "\n";
    ss << "dimension            : " << dimension << "\n";
    ss << "time step            : " << dt_ << "\n";
    ss << "number of time steps : " << timeStepNbr_ << "\n";
    ss << core::to_str(state_->hybridModel_->state);
    ss << "\n";
    ss << "SolverPPC: " << typeid(SolverPPC).name() << "\n";
    return ss.str();
}

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         bool offload>
void Simulator<_dimension, _interp_order, _nbRefinedPart, offload>::initialize()
{
    try
    {
        if (isInitialized)
            std::runtime_error("cannot initialize  - simulator already isInitialized");

        auto& integrator_ = state_->integrator_;

        if (integrator_ == nullptr)
            throw std::runtime_error("Error - Simulator has no integrator");

        integrator_->initialize();

        // if constexpr (!WorkingState_t::is_cpu)
        //     solver::fill_state_from(hostState_->hybridModel_, state_->hybridModel_);
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

    if (core::mpi::any(core::Errors::instance().any()))
    {
        this->dMan.release(); // closes/flushes hdf5 files
        throw std::runtime_error("forcing error");
    }

    isInitialized = true;
}




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         bool offload>
double Simulator<_dimension, _interp_order, _nbRefinedPart, offload>::advance(double dt)
{
    double dt_new     = 0;
    auto& integrator_ = state_->integrator_;

    if (!integrator_)
        throw std::runtime_error("Error - no valid integrator in the simulator");

    try
    {
        PHARE_LOG_SCOPE("Simulator::advance");
        dt_new       = integrator_->advance(dt);
        currentTime_ = ((*timeStamper) += dt);
    }
    catch (std::runtime_error const& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }

    if (core::mpi::any(core::Errors::instance().any()))
    {
        this->dMan.release(); // closes/flushes hdf5 files
        throw std::runtime_error("forcing error");
    }

    return dt_new;
}




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         bool offload>
auto Simulator<_dimension, _interp_order, _nbRefinedPart, offload>::find_model(std::string name)
{
    return std::find(std::begin(modelNames_), std::end(modelNames_), name) != std::end(modelNames_);
}



struct SimulatorMaker
{
    SimulatorMaker(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
        : hierarchy_{hierarchy}
    {
    }

    std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy_;

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    std::unique_ptr<ISimulator> operator()(std::size_t userDim, std::size_t userInterpOrder,
                                           std::size_t userNbRefinedPart, Dimension dimension,
                                           InterpOrder interp_order, NbRefinedPart nbRefinedPart)
    {
        if (userDim == dimension() and userInterpOrder == interp_order()
            and userNbRefinedPart == nbRefinedPart())
        {
            std::size_t constexpr d  = dimension();
            std::size_t constexpr io = interp_order();
            std::size_t constexpr nb = nbRefinedPart();

            auto& theDict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();

#if defined(PHARE_WITH_GPU)
            return std::make_unique<Simulator<d, io, nb, 1>>(theDict, hierarchy_);
#else
            return std::make_unique<Simulator<d, io, nb, 0>>(theDict, hierarchy_);
#endif
        }

        return nullptr;
    }
};


std::unique_ptr<PHARE::ISimulator> getSimulator(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy);


template<typename Simulator>
std::unique_ptr<Simulator> makeSimulator(std::shared_ptr<amr::Hierarchy> const& hierarchy)
{
    return std::make_unique<Simulator>(initializer::PHAREDictHandler::INSTANCE().dict(), hierarchy);
}


} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_H*/
