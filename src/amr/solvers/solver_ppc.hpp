

#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include <thread>

#include <SAMRAI/hier/Patch.h>

#include "initializer/data_provider.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/resources_manager/resources_manager.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/solvers/solver.hpp"

#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/ohm/ohm.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/grid/gridlayout_utils.hpp"


#include "core/numerics/ion_updater/ion_range.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "core/utilities/thread_pool.hpp"

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


    std::size_t static operating_n_threads(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict.contains("threads"))
            return dict["threads"].template to<std::size_t>();
        return 1;
    }

    PHARE::initializer::PHAREDict updaterDict;

public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    explicit SolverPPC(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"PPC"}
        , ohm_{dict["ohm"]}
        , updaterDict{dict["ion_updater"]}
        , n_threads{operating_n_threads(dict)}
    {
        thread_pool_ = std::make_unique<thread_pool>(n_threads);
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

    using IonPopView   = core::IonPopulationView<ParticleArray, VecFieldT, GridLayout>;
    using IonUpdater_t = core::IonUpdater<typename Electromag::view_t, ParticleArray, GridLayout>;
    using RangeSynchrotron_t = PHARE::core::RangeSynchrotron<ParticleArray>;

    using PatchView = std::tuple<GridLayout, typename Electromag::view_t,
                                 std::vector<std::shared_ptr<IonPopView>>>;
    std::vector<PatchView> ion_patch_views;

    std::unique_ptr<thread_pool> thread_pool_;
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

    auto reset_moments = [&]() {
        auto& ions = hybridState.ions;
        for (auto& patch : *level)
        {
            auto _ = resourcesManager.setOnPatch(*patch, ions);
            resetMoments(ions);
        }
    };

    ion_patch_views.clear();
    for (auto& patch : *level)
    {
        auto _      = resourcesManager.setOnPatch(*patch, hybridState);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        ion_patch_views.emplace_back(layout, hybridState.electromag.view(),
                                     IonPopView::make_shared(hybridState.ions));
    }

    predictor1_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);

    saveState_(*level, hybridState.ions, resourcesManager);
    reset_moments();
    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::domain_only);

    predictor2_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);

    restoreState_(*level, hybridState.ions, resourcesManager);
    reset_moments();
    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::all);

    corrector_(*level, hybridModel, fromCoarser, currentTime, newTime);
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
        auto units = PHARE::core::updater_ranges_per_thread(ion_patch_views, n_threads);
        core::abort_if(units.size() != n_threads);

        auto synchrotron
            = std::make_shared<RangeSynchrotron_t>(static_cast<std::uint16_t>(n_threads));

        auto thread_fn = [&](std::uint16_t thread_idx) {
            IonUpdater_t{updaterDict, thread_idx, synchrotron}.updatePopulations(units[thread_idx],
                                                                                 dt, mode);
        };
        for (std::size_t i = 1; i < n_threads; ++i)
            thread_pool_->submit([&, i]() { thread_fn(i); });
        thread_fn(0);
        thread_pool_->wait_for_tasks();
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
        IonUpdater_t::updateIons(ions, layout);
        // no need to update time, since it has been done before
    }
}
} // namespace PHARE::solver


#endif
