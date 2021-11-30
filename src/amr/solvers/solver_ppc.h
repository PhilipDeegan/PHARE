

#ifndef PHARE_SOLVER_PPC_H
#define PHARE_SOLVER_PPC_H

#include <iomanip>

#include <SAMRAI/hier/Patch.h>

#include "initializer/data_provider.h"

#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/hybrid_messenger_info.h"
#include "amr/resources_manager/amr_utils.h"
#include "amr/resources_manager/resources_manager.h"
#include "amr/resources_manager/amr_utils.h"

#include "amr/solvers/solver.h"

#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/boundary_condition/boundary_condition.h"

#include "core/numerics/ampere/ampere.h"
#include "core/numerics/faraday/faraday.h"
#include "core/numerics/ohm/ohm.h"

#include "core/data/particles/particle_array.h"
#include "core/data/vecfield/vecfield.h"
#include "core/data/grid/gridlayout_utils.h"


#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/numerics/ion_updater/ion_range_updater.h"



namespace PHARE::solver
{
// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
class SolverPPC : public ISolver<AMR_Types>
{
private:
    static constexpr auto dimension    = HybridModel::dimension;
    static constexpr auto interp_order = HybridModel::gridlayout_type::interp_order;

    using Electromag       = typename HybridModel::electromag_type;
    using Ions             = typename HybridModel::ions_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename HybridModel::vecfield_type;
    using GridLayout       = typename HybridModel::gridlayout_type;
    using ResourcesManager = typename HybridModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using HybridMessenger  = amr::HybridMessenger<HybridModel>;


    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};


    PHARE::core::Faraday<GridLayout> faraday_;
    PHARE::core::Ampere<GridLayout> ampere_;
    PHARE::core::Ohm<GridLayout> ohm_;
    PHARE::core::IonUpdater<Ions, Electromag, GridLayout> ionUpdater_;

    std::size_t static operating_n_threads(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict.contains("threads"))
            return dict["threads"].template to<std::size_t>();
        return 1;
    }


public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    explicit SolverPPC(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"PPC"}
        , ohm_{dict["ohm"]}
        , ionUpdater_{dict["ion_updater"]}
        , n_threads{operating_n_threads(dict)}
    {
    }

    virtual ~SolverPPC() = default;


    virtual std::string modelName() const override { return HybridModel::model_name; }

    virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    virtual void registerResources(IPhysicalModel_t& model) override;


    virtual void allocate(IPhysicalModel_t& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override;



    virtual void advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                              IPhysicalModel_t& model, IMessenger& fromCoarserMessenger,
                              double const currentTime, double const newTime) override;



private:
    using Messenger = amr::HybridMessenger<HybridModel>;


    void predictor1_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void predictor2_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void corrector_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                    double const currentTime, double const newTime);


    void average_(level_t& level, HybridModel& model);


    void moveIons_(level_t& level, Ions& ions, Electromag& electromag, ResourcesManager& rm,
                   Messenger& fromCoarser, double const currentTime, double const newTime,
                   core::UpdaterMode mode);


    void saveState_(level_t& level, Ions& ions, ResourcesManager& rm);

    void restoreState_(level_t& level, Ions& ions, ResourcesManager& rm);

    /*
    template<typename HybridMessenger>
    void syncLevel(HybridMessenger& toCoarser)
    {
        toCoarser.syncMagnetic(model_.electromag.B);
        toCoarser.syncElectric(model_.electromag.E);
    }*/

    std::size_t n_threads = 1;

    // extend lifespan
    std::unordered_map<std::string, ParticleArray> tmpDomain;
    std::unordered_map<std::string, ParticleArray> patchGhost;

    using StateView_t       = typename HybridModel::StateView_t;
    using IonRangeUpdater_t = PHARE::core::IonRangeUpdater<HybridModel, GridLayout>;
    std::vector<std::shared_ptr<StateView_t>> state_views;

    using RangeSynchrotron_t = PHARE::core::RangeSynchrotron<ParticleArray>;


}; // end solverPPC



// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::registerResources(IPhysicalModel_t& model)
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->registerResources(electromagPred_);
    hmodel.resourcesManager->registerResources(electromagAvg_);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::allocate(IPhysicalModel_t& model,
                                                 SAMRAI::hier::Patch& patch,
                                                 double const allocateTime) const
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
    hmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& modelInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    auto const& Epred = electromagPred_.E;
    auto const& Bpred = electromagPred_.B;

    modelInfo.ghostElectric.emplace_back(Epred);
    modelInfo.ghostMagnetic.emplace_back(Bpred);
    modelInfo.initMagnetic.emplace_back(Bpred);
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::saveState_(level_t& level, Ions& ions, ResourcesManager& rm)
{
    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            tmpDomain[ss.str() + "_" + pop.name()]  = pop.domainParticles();
            patchGhost[ss.str() + "_" + pop.name()] = pop.patchGhostParticles();
        }
    }
}

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::restoreState_(level_t& level, Ions& ions,
                                                      ResourcesManager& rm)
{
    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            pop.domainParticles()     = std::move(tmpDomain[ss.str() + "_" + pop.name()]);
            pop.patchGhostParticles() = std::move(patchGhost[ss.str() + "_" + pop.name()]);
        }
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                     int const levelNumber, IPhysicalModel_t& model,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::advanceLevel");

    auto& hybridModel      = dynamic_cast<HybridModel&>(model);
    auto& hybridState      = hybridModel.state;
    auto& fromCoarser      = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *hybridModel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);

    auto set_state_views = [&]() {
        state_views.clear();
        state_views.reserve(level->getNumberOfPatches());
        for (auto& patch : *level)
        {
            auto _      = resourcesManager.setOnPatch(*patch, hybridState);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            state_views.emplace_back(std::make_shared<StateView_t>(hybridState, layout));
        }
    };

    auto reset_moments = [&]() {
        auto& ions = hybridState.ions;
        for (auto& patch : *level)
        {
            auto _ = resourcesManager.setOnPatch(*patch, ions);
            resetMoments(ions);
        }
    };

    predictor1_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);

    reset_moments();
    saveState_(*level, hybridState.ions, resourcesManager);
    set_state_views();
    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::domain_only);

    predictor2_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);

    reset_moments();
    restoreState_(*level, hybridState.ions, resourcesManager);
    set_state_views();
    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::all);

    corrector_(*level, hybridModel, fromCoarser, currentTime, newTime);


    // return newTime;
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::predictor1_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = hybridState.electromag;
    auto levelNumber       = level.getLevelNumber();


    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.faraday");

        auto& Bpred = electromagPred_.B;
        auto& B     = electromag.B;
        auto& E     = electromag.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, B, E);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, E, Bpred, dt);


            resourcesManager->setTime(Bpred, *patch, newTime);
        }

        fromCoarser.fillMagneticGhosts(Bpred, levelNumber, newTime);
    }



    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.ampere");

        auto& Bpred = electromagPred_.B;
        auto& J     = hybridState.J;


        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, ampere_);
            ampere_(Bpred, J);

            resourcesManager->setTime(J, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(J, levelNumber, newTime);
    }



    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.ohm");

        auto& electrons = hybridState.electrons;
        auto& Bpred     = electromagPred_.B;
        auto& Epred     = electromagPred_.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, Epred, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, Bpred, J, Epred);
            resourcesManager->setTime(Epred, *patch, newTime);
        }

        fromCoarser.fillElectricGhosts(Epred, levelNumber, newTime);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::predictor2_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto levelNumber       = level.getLevelNumber();



    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.faraday");

        auto& Bpred = electromagPred_.B;
        auto& B     = hybridState.electromag.B;
        auto& Eavg  = electromagAvg_.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, B, Eavg);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Eavg, Bpred, dt);

            resourcesManager->setTime(Bpred, *patch, newTime);
        }

        fromCoarser.fillMagneticGhosts(Bpred, levelNumber, newTime);
    }


    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.ampere");

        auto& Bpred = electromagPred_.B;
        auto& J     = hybridState.J;


        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, ampere_);
            ampere_(Bpred, J);

            resourcesManager->setTime(J, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(J, levelNumber, newTime);
    }


    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.ohm");

        auto& electrons = hybridState.electrons;
        auto& Bpred     = electromagPred_.B;
        auto& Epred     = electromagPred_.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, Epred, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, Bpred, J, Epred);
            resourcesManager->setTime(Epred, *patch, newTime);
        }

        fromCoarser.fillElectricGhosts(Epred, levelNumber, newTime);
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModel& model,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::corrector_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto levelNumber       = level.getLevelNumber();

    {
        PHARE_LOG_SCOPE("SolverPPC::corrector_.faraday");

        auto& B    = hybridState.electromag.B;
        auto& Eavg = electromagAvg_.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B, Eavg);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Eavg, B, dt);

            resourcesManager->setTime(B, *patch, newTime);
        }

        fromCoarser.fillMagneticGhosts(B, levelNumber, newTime);
    }



    {
        PHARE_LOG_SCOPE("SolverPPC::corrector_.ohm");

        auto& electrons = hybridState.electrons;
        auto& B         = hybridState.electromag.B;
        auto& E         = hybridState.electromag.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, B, E, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, B, J, E);
            resourcesManager->setTime(E, *patch, newTime);
        }

        fromCoarser.fillElectricGhosts(E, levelNumber, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModel& model)
{
    PHARE_LOG_SCOPE("SolverPPC::average_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& Epred = electromagPred_.E;
    auto& Bpred = electromagPred_.B;
    auto& Bavg  = electromagAvg_.B;
    auto& Eavg  = electromagAvg_.E;
    auto& B     = hybridState.electromag.B;
    auto& E     = hybridState.electromag.E;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, electromagAvg_, electromagPred_,
                                              hybridState.electromag);
        PHARE::core::average(B, Bpred, Bavg);
        PHARE::core::average(E, Epred, Eavg);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, Ions& ions,
                                                  Electromag& electromag, ResourcesManager& rm,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE("SolverPPC::moveIons_");

    std::size_t nbrDomainParticles        = 0;
    std::size_t nbrPatchGhostParticles    = 0;
    std::size_t nbrLevelGhostNewParticles = 0;
    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0;
    for (auto& patch : level)
    {
        auto _ = rm.setOnPatch(*patch, ions);

        for (auto& pop : ions)
        {
            nbrDomainParticles += pop.domainParticles().size();
            nbrPatchGhostParticles += pop.patchGhostParticles().size();
            nbrLevelGhostNewParticles += pop.levelGhostParticlesNew().size();
            nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
            nbrLevelGhostParticles += pop.levelGhostParticles().size();
            nbrPatchGhostParticles += pop.patchGhostParticles().size();

            if (nbrLevelGhostOldParticles < nbrLevelGhostParticles
                and nbrLevelGhostOldParticles > 0)
                throw std::runtime_error("Error - there are less old level ghost particles ("
                                         + std::to_string(nbrLevelGhostOldParticles)
                                         + ") than pushable ("
                                         + std::to_string(nbrLevelGhostParticles) + ")");
        }
    }

    core::abort_if(n_threads == 0);
    auto dt = newTime - currentTime;

    { // syncs on destruct
        auto units = PHARE::core::solver_update_dao_per_thread(state_views, n_threads);
        core::abort_if(units.size() != n_threads);

        RangeSynchrotron_t synchrotron{static_cast<std::uint16_t>(n_threads)};

        auto thread_fn = [&](std::uint16_t idx) {
            for (auto& pop : units[idx])
                IonRangeUpdater_t{synchrotron, idx, "modified_boris"}.updatePopulations(pop, dt,
                                                                                        mode);
        };
        auto threads = PHARE::core::generate(
            [&](auto i) { return std::thread{[&, i]() { thread_fn(i); }}; }, 1, n_threads);
        thread_fn(0);
        for (auto& thread : threads)
            if (thread.joinable())
                thread.join();
    }

    for (auto& patch : level)
    {
        auto _ = rm.setOnPatch(*patch, electromag, ions);
        // this needs to be done before calling the messenger
        rm.setTime(ions, *patch, newTime);
    }

    fromCoarser.fillIonGhostParticles(ions, level, newTime);
    fromCoarser.fillIonMomentGhosts(ions, level, currentTime, newTime);

    for (auto& patch : level)
    {
        auto _      = rm.setOnPatch(*patch, electromag, ions);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        ionUpdater_.updateIons(ions, layout);
        // no need to update time, since it has been done before
    }
}
} // namespace PHARE::solver


#endif
