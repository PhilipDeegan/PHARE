#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------
#include "core/def/phare_mpi.hpp"
#include "core/data/grid/grid_tiles.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"


#include <SAMRAI/hier/TimeInterpolateOperator.h>


namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename GridLayoutT, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
    static std::size_t constexpr dim = GridLayoutT::dimension;
    static_assert(dim > 0 && dim <= 3);

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;
    using FieldGeometry_t  = FieldGeometry<GridLayoutT, PhysicalQuantity>;
    bool static constexpr withGhost{true};

    auto static toFieldBox(auto&&... args) { return FieldGeometry_t::toFieldBox(args...); }

public:
    using GridLayoutImpl = typename GridLayoutT::implT;

    FieldLinearTimeInterpolate()
        : SAMRAI::hier::TimeInterpolateOperator{"FieldLinearTimeInterpolate"}
    {
    }


    virtual ~FieldLinearTimeInterpolate() = default;

    void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& where,
                         SAMRAI::hier::BoxOverlap const& /*overlap*/,
                         SAMRAI::hier::PatchData const& srcDataOld,
                         SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        auto& fieldDataDest = dynamic_cast<FieldDataT&>(destData);

        auto const& fieldDataSrcOld = dynamic_cast<FieldDataT const&>(srcDataOld);
        auto const& fieldDataSrcNew = dynamic_cast<FieldDataT const&>(srcDataNew);

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;
        auto& fieldDest         = fieldDataDest.field;

        double const interpTime = fieldDataDest.getTime();
        double const oldTime    = fieldDataSrcOld.getTime();
        double const newTime    = fieldDataSrcNew.getTime();
        double const alpha      = (interpTime - oldTime) / (newTime - oldTime);
        auto const qty          = fieldDest.physicalQuantity();

        auto const ghostBox
            = toFieldBox(fieldDataDest.getBox(), qty, fieldDataDest.gridLayout, withGhost);

        auto const srcGhostBox
            = toFieldBox(fieldDataSrcNew.getBox(), qty, fieldDataSrcNew.gridLayout, withGhost);

        auto const interpolateBox = toFieldBox(
            where, qty, FieldGeometry_t::layoutFromBox(where, fieldDataDest.gridLayout),
            !withGhost);

        if constexpr (core::is_field_tile_set_v<FieldT>)
        {
            auto const sgbox = phare_box_from<dim>(srcGhostBox);
            auto const dgbox = phare_box_from<dim>(ghostBox);

            for (std::size_t tidx = 0; tidx < fieldSrcOld().size(); ++tidx)
            {
                auto& srcTileOld   = fieldSrcOld()[tidx];
                auto& srcTileNew   = fieldSrcNew()[tidx];
                auto const src_box = srcTileOld.ghost_box();

                if (auto const src_overlap = src_box * sgbox)
                    for (auto& dst_tile : fieldDest())
                    {
                        auto const dst_box = dst_tile.ghost_box();

                        if (auto const dst_overlap = dst_box * dgbox)
                        {
                            auto const finalBox = interpolateBox * samrai_box_from(*dst_overlap);

                            timeInterpolateField(dst_tile(), srcTileOld(), srcTileNew(), qty, alpha,
                                                 samrai_box_from(dst_box), samrai_box_from(src_box),
                                                 finalBox);
                        }
                    }
            }
        }
        else
        {
            auto const finalBox = interpolateBox * ghostBox;
            timeInterpolateField(fieldDest, fieldSrcOld, fieldSrcNew, qty, alpha, ghostBox,
                                 srcGhostBox, finalBox);
        }
    }

    void timeInterpolateField(auto& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                              auto const qty, auto const alpha, auto const ghostBox,
                              auto const srcGhostBox, auto const finalBox) const
    {
        auto const localDestBox = AMRToLocal(finalBox, ghostBox);
        auto const localSrcBox  = AMRToLocal(finalBox, srcGhostBox);

        if constexpr (dim == 1)
        {
            auto const iDestStartX = localDestBox.lower(dirX);
            auto const iDestEndX   = localDestBox.upper(dirX);

            auto const iSrcStartX = localSrcBox.lower(dirX);

            for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
            {
                fieldDest(ix) = (1. - alpha) * fieldSrcOld(ixSrc) + alpha * fieldSrcNew(ixSrc);
            }
        }
        else if constexpr (dim == 2)
        {
            auto const iDestStartX = localDestBox.lower(dirX);
            auto const iDestEndX   = localDestBox.upper(dirX);
            auto const iDestStartY = localDestBox.lower(dirY);
            auto const iDestEndY   = localDestBox.upper(dirY);

            auto const iSrcStartX = localSrcBox.lower(dirX);
            auto const iSrcStartY = localSrcBox.lower(dirY);

            for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
            {
                for (auto iy = iDestStartY, iySrc = iSrcStartY; iy <= iDestEndY; ++iy, ++iySrc)
                {
                    fieldDest(ix, iy) = (1. - alpha) * fieldSrcOld(ixSrc, iySrc)
                                        + alpha * fieldSrcNew(ixSrc, iySrc);
                }
            }
        }
        else if constexpr (dim == 3)
        {
            auto const iDestStartX = localDestBox.lower(dirX);
            auto const iDestEndX   = localDestBox.upper(dirX);
            auto const iDestStartY = localDestBox.lower(dirY);
            auto const iDestEndY   = localDestBox.upper(dirY);
            auto const iDestStartZ = localDestBox.lower(dirZ);
            auto const iDestEndZ   = localDestBox.upper(dirZ);

            auto const iSrcStartX = localSrcBox.lower(dirX);
            auto const iSrcStartY = localSrcBox.lower(dirY);
            auto const iSrcStartZ = localSrcBox.lower(dirZ);

            for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
            {
                for (auto iy = iDestStartY, iySrc = iSrcStartY; iy <= iDestEndY; ++iy, ++iySrc)
                {
                    for (auto iz = iDestStartZ, izSrc = iSrcStartZ; iz <= iDestEndZ; ++iz, ++izSrc)
                    {
                        fieldDest(ix, iy, iz) = (1. - alpha) * fieldSrcOld(ixSrc, iySrc, izSrc)
                                                + alpha * fieldSrcNew(ixSrc, iySrc, izSrc);
                    }
                }
            }
        }
    }
};

} // namespace PHARE::amr

#endif
