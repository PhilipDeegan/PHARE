#ifndef PHARE_FIELD_DATA_COARSEN_HPP
#define PHARE_FIELD_DATA_COARSEN_HPP


#include "core/def/phare_mpi.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"

#include "default_field_coarsener.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/IntVector.h>


namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    //
    template<typename GridLayoutT, typename FieldT, typename FieldCoarsenerPolicy,
             typename PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity())>
    /**
     * @brief The FieldCoarsenOperator class
     */
    class FieldCoarsenOperator : public SAMRAI::hier::CoarsenOperator
    {
        static constexpr std::size_t n_ghosts
            = GridLayoutT::template nbrGhosts<core::QtyCentering, core::QtyCentering::dual>();

    public:
        static constexpr std::size_t dimension = GridLayoutT::dimension;
        using FieldDataT                       = FieldData<GridLayoutT, FieldT>;

        FieldCoarsenOperator()
            : SAMRAI::hier::CoarsenOperator("FieldDataCoarsenOperator")
        {
        }

        FieldCoarsenOperator(FieldCoarsenOperator const&)            = delete;
        FieldCoarsenOperator(FieldCoarsenOperator&&)                 = delete;
        FieldCoarsenOperator& operator=(FieldCoarsenOperator const&) = delete;
        FieldCoarsenOperator& operator=(FieldCoarsenOperator&&)      = delete;


        virtual ~FieldCoarsenOperator() = default;




        /** @brief return the priority of the operator
         *  this return 0, meaning that this operator
         * have the most priority
         */
        int getOperatorPriority() const override { return 0; }




        /** @brief Return the stencil width associated with the coarsening operator.
         *
         *  The SAMRAI transfer routines guarantee that the source patch will contain
         * sufficient ghostCell data surrounding the interior to satisfy the stencil
         * width requirements for each coarsening operator.
         *
         * In our case, we allow a RF up to 10, so having 5 ghost width is sufficient
         *
         */
        SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
        {
            return SAMRAI::hier::IntVector{dim, 2};
        }




        /** @brief given a coarseBox, coarse data from the fine patch on the intersection of this
         * box and the box of the destination (the box of the coarse patch).
         *
         * This method will extract fieldData from the two patches, and then
         * get the Field and GridLayout encapsulated into the fieldData.
         * With the help of FieldGeometry, transform the coarseBox to the correct index.
         * After that we can now create FieldCoarsen with the indexAndWeight implementation
         * selected. Finnaly loop over the indexes in the box, and apply the coarsening defined in
         * FieldCoarsen operator
         *
         */
        void coarsen(SAMRAI::hier::Patch& destinationPatch, SAMRAI::hier::Patch const& sourcePatch,
                     int const destinationId, int const sourceId,
                     SAMRAI::hier::Box const& coarseBox,
                     SAMRAI::hier::IntVector const& ratio) const override
        {
            auto& destinationField        = FieldDataT::getField(destinationPatch, destinationId);
            auto const& sourceField       = FieldDataT::getField(sourcePatch, sourceId);
            auto const& sourceLayout      = FieldDataT::getLayout(sourcePatch, sourceId);
            auto const& destinationLayout = FieldDataT::getLayout(destinationPatch, destinationId);

            // we assume that quantity are the same
            // note that an assertion will be raised
            // in coarseIt operator
            auto const& qty = destinationField.physicalQuantity();

            bool const withGhost{true};

            // We get different boxes : destination , source, restrictBoxes
            // and transform them in the correct indexing.
            auto const destinationBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                destinationPatch.getBox(), qty, destinationLayout, withGhost);

            auto const sourceBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                sourcePatch.getBox(), qty, sourceLayout, withGhost);

            auto const coarseLayout = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(
                coarseBox, destinationLayout);

            auto const coarseFieldBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                coarseBox, qty, coarseLayout, !withGhost);

            // finnaly we compute the intersection
            auto const intersectionBox = destinationBox * coarseFieldBox;
            auto const box             = phare_box_from<dimension>(intersectionBox);

            if constexpr (core::is_field_tile_set_v<FieldT>)
            {
                auto const phare_box = phare_box_from<dimension>(intersectionBox);

                for (auto& dst_tile : destinationField())
                {
                    auto const dst_box = dst_tile.field_box();
                    if (auto const dst_overlap = dst_box * box)
                    {
                        for (auto const& src_tile : sourceField())
                        {
                            auto const src_box = src_tile.field_box();
                            if (auto const src_overlap = src_box * refine_box(box))
                            {
                                if (auto const& overlap = *dst_overlap * coarsen_box(*src_overlap))
                                {
                                    FieldCoarsenerPolicy coarsener{
                                        destinationLayout.centering(qty),
                                        samrai_box_from(src_tile.ghost_box()),
                                        samrai_box_from(dst_tile.ghost_box()), ratio};

                                    for (auto const& bix : *overlap)
                                        coarsener(src_tile(), dst_tile(), *bix);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // We can now create the coarsening operator
                FieldCoarsenerPolicy coarsener{destinationLayout.centering(qty), sourceBox,
                                               destinationBox, ratio};
                for (auto const& bix : box)
                    coarsener(sourceField, destinationField, *bix);
            }
        }
    };
} // namespace amr
} // namespace PHARE


#endif
