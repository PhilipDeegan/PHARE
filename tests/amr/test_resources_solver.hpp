#ifndef PHARE_TEST_AMR_RESOURCES_SOLVER_HPP
#define PHARE_TEST_AMR_RESOURCES_SOLVER_HPP

#include "core/models/quantities/hybrid_quantities.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <memory>
#include <stdexcept>


namespace PHARE::amr::test
{

// stands in for SolverPPC in tests that only need the hierarchy set up and filled by messengers:
// registers, allocates and announces to the messenger the same quantities, with the same names,
// as SolverPPC, but never advances. Using SolverPPC instead emits its vtable, which instantiates
// the whole step (ion updater, Ohm, Faraday, Ampere, reflux) per permutation, unused, and can
// get test binaries killed for memory when linking
template<typename HybridModel, typename AMR_Types>
class ResourcesSolver : public solver::ISolver<AMR_Types>
{
    using Super            = solver::ISolver<AMR_Types>;
    using Electromag       = HybridModel::electromag_type;
    using VecFieldT        = HybridModel::vecfield_type;
    using IPhysicalModel_t = Super::IPhysicalModel_t;
    using IMessenger_t     = amr::IMessenger<IPhysicalModel_t>;

public:
    using patch_t     = Super::patch_t;
    using level_t     = Super::level_t;
    using hierarchy_t = Super::hierarchy_t;

    ResourcesSolver()
        : Super{"PPC"} // same name as SolverPPC, quantity names derive from it
    {
    }

    std::string modelName() const override { return HybridModel::model_name; }

    void registerResources(IPhysicalModel_t& model) override
    {
        auto& rm = *dynamic_cast<HybridModel&>(model).resourcesManager;
        rm.registerResources(electromagPred_);
        rm.registerResources(electromagAvg_);
        rm.registerResources(Bold_);
        rm.registerResources(fluxSumE_);
    }

    void allocate(IPhysicalModel_t& model, patch_t& patch, double const allocateTime) const override
    {
        auto& rm = *dynamic_cast<HybridModel&>(model).resourcesManager;
        rm.allocate(electromagPred_, patch, allocateTime);
        rm.allocate(electromagAvg_, patch, allocateTime);
        rm.allocate(Bold_, patch, allocateTime);
        rm.allocate(fluxSumE_, patch, allocateTime);
    }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override
    {
        auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

        auto const& Eavg  = electromagAvg_.E;
        auto const& Bpred = electromagPred_.B;

        hybridInfo.ghostElectric.emplace_back(Eavg.name());
        hybridInfo.initMagnetic.emplace_back(Bpred.name());
        hybridInfo.ghostMagnetic.emplace_back(Bpred.name());
        hybridInfo.refluxElectric  = Eavg.name();
        hybridInfo.fluxSumElectric = fluxSumE_.name();
    }

    void prepareStep(IPhysicalModel_t&, level_t&, double const) override { unsupported(); }
    void accumulateFluxSum(IPhysicalModel_t&, level_t&, double const) override { unsupported(); }
    void resetFluxSum(IPhysicalModel_t&, level_t&) override { unsupported(); }

    void reflux(IPhysicalModel_t&, level_t&, IMessenger_t&, double const) override
    {
        unsupported();
    }

    void advanceLevel(hierarchy_t const&, int const, IPhysicalModel_t&, IMessenger_t&, double const,
                      double const) override
    {
        unsupported();
    }

    double computeStableDt(IPhysicalModel_t&, level_t&, solver::CFLNumbers const&) override
    {
        unsupported();
        return 0;
    }

private:
    [[noreturn]] static void unsupported()
    {
        throw std::runtime_error("ResourcesSolver does not advance, use SolverPPC");
    }

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    VecFieldT Bold_{this->name() + "_Bold", core::HybridQuantity::Vector::B};
    VecFieldT fluxSumE_{this->name() + "_fluxSumE", core::HybridQuantity::Vector::E};
};

} // namespace PHARE::amr::test

#endif /* PHARE_TEST_AMR_RESOURCES_SOLVER_HPP */
