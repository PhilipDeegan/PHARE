#ifndef PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/utilities/constants.hpp"

#include "amr/data/field/field_data_traits.hpp"
#include "amr/data/tensorfield/tensor_field_data_traits.hpp"

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <array>
#include <cassert>
#include <memory>
#include <span>
#include <stdexcept>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename ResMan, typename ScalarOrTensorFieldDataT, typename BoundaryManagerT>
    requires(IsFieldData<ScalarOrTensorFieldDataT> || IsTensorFieldData<ScalarOrTensorFieldDataT>)
class FieldRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    static constexpr bool is_scalar = IsFieldData<ScalarOrTensorFieldDataT>;
    static constexpr bool is_tensor = IsTensorFieldData<ScalarOrTensorFieldDataT>;

    using FieldGeometryT         = FieldGeometrySelector<ScalarOrTensorFieldDataT, is_scalar>::type;
    using GridLayoutT            = ScalarOrTensorFieldDataT::gridlayout_type;
    using GridT                  = ScalarOrTensorFieldDataT::grid_type;
    using FieldT                 = GridT::field_type;
    using PatchGeometry          = SAMRAI::hier::PatchGeometry;
    using CartesianPatchGeometry = SAMRAI::geom::CartesianPatchGeometry;
    using BoundaryType           = core::Bdry::Type;

    static constexpr std::size_t dimension      = ScalarOrTensorFieldDataT::dimension;
    static constexpr int boundaries_codimension = 1; // only handle codimension 1 boundaries at
                                                     ///< the momement

    FieldRefinePatchStrategy(ResMan& resourcesManager, BoundaryManagerT& boundaryManager)
        : rm_{resourcesManager}
        , boundaryManager_{boundaryManager}
        , data_id_{-1}
    {
    }

    void assertIDsSet() const
    {
        assert(data_id_ >= 0 && "FieldRefinePatchStrategy: IDs must be registered before use");
    }

    void registerIDs(int const field_id) { data_id_ = field_id; }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       SAMRAI::hier::IntVector const& ghost_width_to_fill) override
    {
        // Retrieve the FieldData ...
        // std::shared_ptr<SAMRAI::hier::PatchData> patchData = patch.getPatchData(data_id_);
        // if (patchData == nullptr)
        //     throw std::runtime_error("no patch data for the corresponding id "
        //                              + std::to_string(data_id_) + " on patch "
        //                              + std::to_string(patch.getLocalId().getValue()));
        // std::shared_ptr<ScalarOrTensorFieldDataT> fieldData =
        // std::dynamic_pointer_cast<ScalarOrTensorFieldDataT>(patchData); if (!fieldData)
        //     throw std::runtime_error("cannot cast to FieldData");

        // ... to retrieve the gridLayout.
        GridLayoutT const& gridLayout = ScalarOrTensorFieldDataT::getLayout(patch, data_id_);

        // consistency check on the number of ghosts
        // SAMRAI::hier::IntVector dataGhostWidths = patchData->getGhostCellWidth();
        if (ghost_width_to_fill != gridLayout.nbrGhosts())
            throw std::runtime_error("Error - inconsistent ghost cell widths");

        // no check this is a valid cast
        std::shared_ptr<CartesianPatchGeometry> patchGeom
            = std::static_pointer_cast<CartesianPatchGeometry>(patch.getPatchGeometry());

        std::vector<SAMRAI::hier::BoundaryBox> const& boundaries
            = patchGeom->getCodimensionBoundaries(boundaries_codimension);

        std::span<GridT> grids_view;
        // retrieve Field from FieldData
        if constexpr (is_scalar)
        {
            GridT& grid = ScalarOrTensorFieldDataT::getField(patch, data_id_);
            grids_view  = {std::addressof(grid), 1};
            // FieldT& field = *(&grid);
        }
        else if constexpr (is_tensor)
        {
            constexpr size_t N          = ScalarOrTensorFieldDataT::N;
            std::array<GridT, N>& grids = ScalarOrTensorFieldDataT::getFields(patch, data_id_);
            grids_view                  = grids;
        };


        // must be retrieved to pass as argument to patchGeom->getBoundaryFillBox later
        SAMRAI::hier::Box const& patch_box = patch.getBox();

        for (SAMRAI::hier::BoundaryBox const& bBox : boundaries)
        {
            // Boundary definitions in PHARE matches those of SAMRAI
            core::BdryLoc::Type const bLoc
                = static_cast<core::BdryLoc::Type>(bBox.getLocationIndex());

            // Get ghost cells as a SAMRAI Box, and extend it in case of dual staggering.
            SAMRAI::hier::Box samraiBoxToFill
                = patchGeom->getBoundaryFillBox(bBox, patch_box, ghost_width_to_fill);
            for (auto const& grid : grids_view)
            {
                FieldT field = *(&grid);

                // take staggering into account to compute the box to fill
                SAMRAI::hier::Box const samraiFieldBoxToFill = FieldGeometryT::toFieldBox(
                    samraiBoxToFill, field.physicalQuantity(), gridLayout);

                // turn SAMRAI box into a local PHARE::core::Box
                auto AMRBoxToFill = phare_box_from<dimension>(samraiFieldBoxToFill);
                auto localBox     = gridLayout.AMRToLocal(AMRBoxToFill);

                auto bc = boundaryManager_.getBoundaryCondition(bLoc, field.physicalQuantity());
                bc->apply(field, localBox, gridLayout, fill_time);
            }
        };
    }

    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector(dim, 1); // hard-coded 0th order base interpolation
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
    }


    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
    }


    static auto isNewFineFace(auto const& amrIdx, auto const dir)
    {
        // amr index can be negative so test !=0 and not ==1
        // to see if this is odd or even
        return amrIdx[dir] % 2 != 0;
    }


protected:
    ResMan& rm_;
    BoundaryManagerT& boundaryManager_;
    int data_id_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
