
#ifndef PHARE_SOLVER_INCLUDE_HPP
#define PHARE_SOLVER_INCLUDE_HPP

#include "phare_core.hpp"
#include "phare_amr.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver_ppc.hpp"
#include "amr/level_initializer/level_initializer.hpp"
#include "amr/level_initializer/level_initializer_factory.hpp"
#include "amr/multiphysics_integrator.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"


namespace PHARE::solver
{

template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_, //
         auto L_ = core::LayoutMode::AoSMapped, auto A_ = AllocatorMode::CPU>
struct PHARE_Types
{
    static_assert(std::is_same_v<decltype(L_), core::LayoutMode>);
    static_assert(std::is_same_v<decltype(A_), AllocatorMode>);

    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    // core deps
    using core_types   = core::PHARE_Types<dimension, interp_order, L_, A_>;
    using VecField_t   = typename core_types::VecField_t;
    using Grid_t       = typename core_types::Grid_t;
    using Electromag_t = typename core_types::Electromag_t;
    using Ions_t       = typename core_types::Ions_t;
    using GridLayout_t = typename core_types::GridLayout_t;
    using Electrons_t  = typename core_types::Electrons_t;
    // core deps

    using HybridModel_t
        = HybridModel<GridLayout_t, Electromag_t, Ions_t, Electrons_t, amr::SAMRAI_Types, Grid_t>;
    using MHDModel_t                = MHDModel<GridLayout_t, VecField_t, amr::SAMRAI_Types, Grid_t>;
    using SolverPPC_t               = SolverPPC<HybridModel_t, amr::SAMRAI_Types>;
    using SolverMHD_t               = SolverMHD<MHDModel_t, amr::SAMRAI_Types>;
    using LevelInitializerFactory_t = LevelInitializerFactory<HybridModel_t>;

    // amr deps
    using amr_types        = amr::PHARE_Types<dimension, interp_order, nbRefinedPart, L_, A_>;
    using RefinementParams = typename amr_types::RefinementParams;

    using MessengerFactory // = amr/solver bidirectional dependency
        = amr::MessengerFactory<MHDModel_t, HybridModel_t, RefinementParams>;
    // amr deps

    using MultiPhysicsIntegrator_t
        = MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t, amr::SAMRAI_Types>;
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_INCLUDE_HPP
