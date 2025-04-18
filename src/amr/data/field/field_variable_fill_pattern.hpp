#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP


#include "core/def/phare_mpi.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include <core/hybrid/hybrid_quantities.hpp>

#include <amr/utilities/box/amr_box.hpp>
#include "amr/data/field/refine/field_refine_operator.hpp"

#include <SAMRAI/pdat/CellOverlap.h>
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <cassert>

namespace PHARE::amr
{
/*
  This class is used from multiple schedules
  To know which schedule we are coming from, we have `std::optional<bool> opt_overwrite_interior_`

  the modes are :

    1. To synchronize primal nodes on shared patch borders
        e.g. hybrid_hybrid_messenger_strategy.hpp
        HybridHybridMessengerStrategy::magneticSharedNodes_

        in this case, the fillPattern is constructed
        with "std::optional<bool> opt_overwrite_interior_ == std::nullopt",
        we set the forwarding flag of "bool overwrite_interior" to true by default
        and it is only set to false for one of the 2 patches involved in the overlap
        so that only one process assigns its value to the shared border node
        We also remove the exclusive interior of the src patch to isolate only shared primal
  nodes.

    2. To synchronize pure ghost values from src domain values
    in that case, the optional is set to "false" and overwrite_interior takes this value
    none of the two patches overwrites the shared border nodes and only pure ghost nodes are
    filled.

  Notes on shared-node overwrite interior: https://github.com/LLNL/SAMRAI/issues/170

*/
// This class is mostly a copy of BoxGeometryVariableFillPattern
template<std::size_t dimension>
class FieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    constexpr static std::size_t dim = dimension;

public:
    FieldFillPattern(std::optional<bool> overwrite_interior)
        : opt_overwrite_interior_{overwrite_interior}
    {
    }

    static auto make_shared(std::shared_ptr<SAMRAI::hier::RefineOperator> const& samrai_op)
    {
        auto const& op = dynamic_cast<AFieldRefineOperator const&>(*samrai_op);

        if (op.node_only)
            return std::make_shared<FieldFillPattern<dim>>(std::nullopt);

        return std::make_shared<FieldFillPattern<dim>>(false);
    }


    virtual ~FieldFillPattern() {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const fn_overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
        NULL_USE(dst_patch_box);
#endif
        TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

        bool overwrite_interior = true; // replace func param
        assert(fn_overwrite_interior == overwrite_interior);

        if (opt_overwrite_interior_) // not node only
        {
            // this sets overwrite_interior to false
            overwrite_interior = *opt_overwrite_interior_;
        }

        // opt_overwrite_interior_ is nullopt : assume primal node shared border schedule
        else
        {
            // cast into the Base class to get the pureInteriorFieldBox method
            // see field_geometry.hpp for more explanations about why this base class exists
            auto& dst_cast = dynamic_cast<FieldGeometryBase<dimension> const&>(dst_geometry);
            auto& src_cast = dynamic_cast<FieldGeometryBase<dimension> const&>(src_geometry);

            if (src_cast.patchBox.getGlobalId().getOwnerRank()
                != dst_cast.patchBox.getGlobalId().getOwnerRank())
                overwrite_interior
                    = src_cast.patchBox.getGlobalId() > dst_cast.patchBox.getGlobalId();

            auto basic_overlap = dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
                                                               overwrite_interior, transformation);
            auto& overlap      = dynamic_cast<FieldOverlap const&>(*basic_overlap);

            auto destinationBoxes = overlap.getDestinationBoxContainer();
            destinationBoxes.removeIntersections(src_cast.pureInteriorFieldBox());

            return std::make_shared<FieldOverlap>(destinationBoxes, overlap.getTransformation());
        }

        // overwrite_interior is always false here
        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
                                             transformation);
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    FieldFillPattern(FieldFillPattern const&)            = delete;
    FieldFillPattern& operator=(FieldFillPattern const&) = delete;

    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        TBOX_ERROR("getStencilWidth() should not be\n"
                   << "called.  This pattern creates overlaps based on\n"
                   << "the BoxGeometry objects and is not restricted to a\n"
                   << "specific stencil.\n");

        /*
         * Dummy return value that will never get reached.
         */
        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(1));
    }

    /*
     *************************************************************************
     *
     * Compute BoxOverlap that specifies data to be filled by refinement
     * operator.
     *
     *************************************************************************
     */
    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        NULL_USE(node_fill_boxes);

        /*
         * For this (default) case, the overlap is simply the intersection of
         * fill_boxes and data_box.
         */
        SAMRAI::hier::Transformation transformation(
            SAMRAI::hier::IntVector::getZero(patch_box.getDim()));

        SAMRAI::hier::BoxContainer overlap_boxes(fill_boxes);
        overlap_boxes.intersectBoxes(data_box);
        return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);
    }

    std::optional<bool> opt_overwrite_interior_{nullptr};
};




template<typename Gridlayout_t> // ASSUMED ALL PRIMAL!
class FieldGhostInterpOverlapFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    std::size_t constexpr static dim         = Gridlayout_t::dimension;
    bool constexpr static overwrite_interior = false;

public:
    FieldGhostInterpOverlapFillPattern() {}
    ~FieldGhostInterpOverlapFillPattern() override {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const fn_overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        PHARE_LOG_SCOPE(3, "FieldGhostInterpOverlapFillPattern::calculateOverlap");

        if (phare_box_from<dim>(dst_patch_box) == phare_box_from<dim>(src_mask))
            return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);

        auto& dst_cast
            = dynamic_cast<FieldGeometry<Gridlayout_t, core::HybridQuantity::Scalar> const&>(
                dst_geometry);

        auto dst_domain    = phare_box_from<dim>(dst_patch_box);
        auto dst_ghost_box = grow(dst_domain, Gridlayout_t::nbrGhosts());
        dst_ghost_box.upper += 1;

        auto bsrc_mask = grow(phare_box_from<dim>(src_mask), Gridlayout_t::nbrGhosts());
        bsrc_mask.upper += 1;

        auto shifted   = core::shift(bsrc_mask, as_point<dim>(transformation));
        auto intersect = *(dst_ghost_box * shifted);

        SAMRAI::hier::BoxContainer dest;
        dest.push_back(samrai_box_from(intersect));
        dest.removeIntersections(samrai_box_from(shrink(
            dst_ghost_box, Gridlayout_t::nbrGhosts() + 1 + (Gridlayout_t::interp_order > 1))));

        return std::make_shared<FieldOverlap>(dest, transformation);
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        throw std::runtime_error("never called");
    }


    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        throw std::runtime_error("no refinement supported or expected");
    }
};



} // namespace PHARE::amr

#endif /* PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_H */
