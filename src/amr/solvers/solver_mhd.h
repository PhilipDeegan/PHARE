
#ifndef PHARE_SOLVER_MHD_H
#define PHARE_SOLVER_MHD_H

#include "amr/messengers/mhd_messenger_info.h"
#include "amr/physical_models/physical_model.h"
#include "amr/solvers/solver.h"

namespace PHARE
{
namespace solver
{
    template<typename MHDModel, typename AMR_Types>
    class SolverMHD : public ISolver<AMR_Types, typename MHDModel::Float>
    {
        using Float            = typename MHDModel::Float;
        using IPhysicalModel_t = IPhysicalModel<AMR_Types, Float>;

    public:
        SolverMHD()
            : ISolver<AMR_Types, Float>{"MHDSolver"}
        {
        }


        virtual ~SolverMHD() = default;

        virtual std::string modelName() const override { return MHDModel::model_name; }


        virtual void
        fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
        {
        }


        virtual void registerResources(IPhysicalModel_t& /*model*/) override {}

        // TODO make this a resourcesUser
        virtual void allocate(IPhysicalModel_t& /*model*/, SAMRAI::hier::Patch& /*patch*/,
                              Float const /*allocateTime*/) const override
        {
        }

        virtual void
        advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                     int const /*levelNumber*/, IPhysicalModel_t& /*model*/,
                     amr::IMessenger<IPhysicalModel_t>& /*fromCoarser*/,
                     const Float /*currentTime*/, const Float /*newTime*/) override
        {
        }
    };
} // namespace solver
} // namespace PHARE


#endif
