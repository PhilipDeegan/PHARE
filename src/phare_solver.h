
#ifndef PHARE_SOLVER_INCLUDE_H
#define PHARE_SOLVER_INCLUDE_H

#include "phare_amr.h"

#include "amr/solvers/solver.h"
#include "amr/solvers/solver_mhd.h"
#include "amr/solvers/solver_ppc.h"

#include "amr/level_initializer/level_initializer.h"
#include "amr/level_initializer/level_initializer_factory.h"
#include "amr/multiphysics_integrator.h"
#include "amr/physical_models/hybrid_model.h"
#include "amr/physical_models/mhd_model.h"
#include "amr/physical_models/physical_model.h"

#if defined(HAVE_RAJA) and defined(HAVE_UMPIRE)
#include "amr/solvers/llnl/solver_ppc.h"
#elif defined(PHARE_WITH_GPU)
#include "amr/solvers/gpu/solver_ppc.h"
#endif

namespace PHARE::solver
{
template<std::size_t dimension, std::size_t interp_order, std::size_t nbRefinedPart>
struct CPU_Types
{
    using core_types = PHARE::core::CPU_Types<dimension, interp_order>;

    using VecField_t   = typename core_types::VecField_t;
    using Electromag_t = typename core_types::Electromag_t;
    using Ions_t       = typename core_types::Ions_t;
    using GridLayout_t = typename core_types::GridLayout_t;
    using Electrons_t  = typename core_types::Electrons_t;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using MHDModel_t = PHARE::solver::MHDModel<GridLayout_t, VecField_t, PHARE::amr::SAMRAI_Types>;
    using HybridModel_t = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t,
                                                     Electrons_t, PHARE::amr::SAMRAI_Types>;


    using SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;


    using LevelInitializerFactory_t = PHARE::solver::LevelInitializerFactory<HybridModel_t>;


    // amr deps
    using amr_types        = PHARE::amr::PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using RefinementParams = typename amr_types::RefinementParams;

    using MessengerFactory // = amr/solver bidirectional dependency
        = PHARE::amr::MessengerFactory<MHDModel_t, HybridModel_t, RefinementParams>;
    // amr deps

    using MultiPhysicsIntegrator
        = PHARE::solver::MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t,
                                                PHARE::amr::SAMRAI_Types>;
};


template<typename HybridModel_t, std::size_t nbRefinedPart, bool offload>
auto constexpr solver_select()
{
    auto constexpr dimension    = HybridModel_t::dimension;
    auto constexpr interp_order = HybridModel_t::gridlayout_type::interp_order;
    using CPU_SolverPPC_t = typename CPU_Types<dimension, interp_order, nbRefinedPart>::SolverPPC_t;

#if defined(HAVE_RAJA) and defined(HAVE_UMPIRE)
    using GPU_SolverPPC_t = PHARE::solver::llnl::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    return static_cast<std::conditional_t<offload, GPU_SolverPPC_t, CPU_SolverPPC_t>*>(nullptr);

#elif defined(PHARE_WITH_GPU)
    using GPU_SolverPPC_t = PHARE::solver::gpu::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    return static_cast<std::conditional_t<offload, GPU_SolverPPC_t, CPU_SolverPPC_t>*>(nullptr);

#elif defined(HAVE_RAJA) or defined(HAVE_UMPIRE)
#error // invalid, both RAJA and UMPIRE are required together.
#else

    return static_cast<CPU_SolverPPC_t*>(nullptr);
#endif
}

template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_,
         bool offload = false>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    // core deps
    using core_types   = PHARE::core::PHARE_Types<dimension, interp_order>;
    using VecField_t   = typename core_types::VecField_t;
    using Electromag_t = typename core_types::Electromag_t;
    using Ions_t       = typename core_types::Ions_t;
    using GridLayout_t = typename core_types::GridLayout_t;
    using Electrons_t  = typename core_types::Electrons_t;
    // core deps

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel_t  = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t,
                                                     Electrons_t, PHARE::amr::SAMRAI_Types>;
    using MHDModel_t = PHARE::solver::MHDModel<GridLayout_t, VecField_t, PHARE::amr::SAMRAI_Types>;

    using SolverPPC_t
        = std::decay_t<decltype(*solver_select<HybridModel_t, nbRefinedPart, offload>())>;

    using SolverMHD_t = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types>;

    // use CPU_Types::LevelInitializerFactory_t - this next one may be for device mem which is
    // copied from cpu mem using
    // LevelInitializerFactory_t = PHARE::solver::LevelInitializerFactory<HybridModel_t>;

    // amr deps
    using amr_types        = PHARE::amr::PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using RefinementParams = typename amr_types::RefinementParams;

    using MessengerFactory // = amr/solver bidirectional dependency
        = PHARE::amr::MessengerFactory<MHDModel_t, HybridModel_t, RefinementParams>;
    // amr deps

    using MultiPhysicsIntegrator = PHARE::solver::MultiPhysicsIntegrator<
        MessengerFactory, PHARE::solver::LevelInitializerFactory<HybridModel_t>,
        PHARE::amr::SAMRAI_Types>;
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_INCLUDE_H
