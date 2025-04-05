#ifndef PHARE_PARTICLES_DATA_SPLIT_HPP
#define PHARE_PARTICLES_DATA_SPLIT_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp"

// #include "core/utilities/constants.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"
#include "core/data/particles/particle_array_exporter.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"


#include "amr/data/particles/particles_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "split.hpp"
#include "amr/amr_constants.hpp"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/CellOverlap.h>

#include <tuple>
// #include <functional>


namespace PHARE
{
namespace amr
{
    enum class ParticlesDataSplitType {
        coarseBoundary,
        interior,
        coarseBoundaryOld,
        coarseBoundaryNew
    };

}
} // namespace PHARE

namespace PHARE::amr
{



template<typename Iterator>
NO_DISCARD auto toFineGrid(Iterator iterator)
{
    constexpr auto dim   = Iterator::dimension;
    constexpr auto ratio = PHARE::amr::refinementRatio;

    core::ParticlePosition<dim> finer{iterator.iCell(), iterator.delta()};
    auto& [iCell, delta] = finer;
    for (size_t iDim = 0; iDim < dim; ++iDim)
    {
        auto const fineDelta    = delta[iDim] * ratio;
        auto const fineDeltaInt = static_cast<int>(fineDelta);
        iCell[iDim]             = iCell[iDim] * ratio + fineDeltaInt;
        delta[iDim]             = fineDelta - fineDeltaInt;
    }

    return finer;
}



template<typename ParticleArray, ParticlesDataSplitType splitType, typename Splitter,
         auto ParticleType_v = ParticleType::Domain>
struct ParticlesRefining
{
    static constexpr auto dim           = Splitter::dimension;
    static constexpr auto interpOrder   = Splitter::interp_order;
    static constexpr auto nbRefinedPart = Splitter::nbRefinedPart;

    ParticlesData<ParticleArray>& srcParticlesData;
    ParticlesData<ParticleArray>& destParticlesData;


    // the particle refine operator's job is to fill either domain (during initialization of
    // new patches) or coarse to fine boundaries (during advance), so we need references to
    // these arrays on the destination. We don't fill ghosts with this operator, they are
    // filled from exchanging with neighbor patches.
    ParticleArray& destParticles = pickDestParticles();

    Splitter split{};

    // the source PatchData is a possible restriction of a "real" patchdata
    // so that it is the closest from the destination boxes
    // if all particles from the original source patchdata are in "domainParticles"
    // they can now be found in either domain of ghost particle arrays of this
    // temporary restriction "source" patchData
    // therefore we need references to the domain and ghost particle arrays
    std::array<ParticleArray const*, 2> particlesArrays{&srcParticlesData.domainParticles,
                                                        &srcParticlesData.patchGhostParticles};

    auto& pickDestParticles()
    {
        bool constexpr putParticlesInCoarseBoundary
            = splitType == ParticlesDataSplitType::coarseBoundary
              || splitType == ParticlesDataSplitType::coarseBoundaryOld
              || splitType == ParticlesDataSplitType::coarseBoundaryNew;

        if constexpr (putParticlesInCoarseBoundary)
        {
            if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
                return destParticlesData.levelGhostParticles;

            else if constexpr (splitType == ParticlesDataSplitType::coarseBoundaryOld)
                return destParticlesData.levelGhostParticlesOld;

            else if constexpr (splitType == ParticlesDataSplitType::coarseBoundaryNew)
                return destParticlesData.levelGhostParticlesNew; /*
            else
                compile error  */
        }

        else
            return destParticlesData.domainParticles;
    }

    void _forBox(core::Box<int, dim> const& destinationBox)
    {
        using ArrayParticleArray = typename ParticleArray::template array_type<nbRefinedPart>;
        auto const splitBox      = getSplitBox(destinationBox);

        auto const per_particle = [&](auto const& particle) {
            auto refined_info = std::tuple<std::uint16_t, ArrayParticleArray>{0, {}};
            auto& [p_count, refinedParticles] = refined_info;
            auto const particleRefinedPos     = toFineGrid(particle);
            split(particleRefinedPos, particle, refinedParticles);
            p_count = ParticleArrayPartitioner<ArrayParticleArray>{refinedParticles}(destinationBox)
                          .size();
            return refined_info;
        };

        auto const refiner = [&](auto const& particle) { return toFineGrid(particle); };

        for (auto const& sourceParticlesArray : particlesArrays)
            export_refined_particles( //
                *sourceParticlesArray, destParticles, splitBox, refiner, per_particle);
    }

    void forBoxes(SAMRAI::hier::BoxContainer const& boxes)
    {
        std::size_t const curr_size = destParticles.size();
        for (auto const& box : boxes)
            _forBox(phare_box_from<dim>(box));
        _finish(curr_size);
    }

    template<typename BoxContainer_t>
    void forBoxes(BoxContainer_t const& boxes)
    {
        std::size_t const curr_size = destParticles.size();
        for (auto const& box : boxes)
            _forBox(box);
        _finish(curr_size);
    }

    void _finish(std::size_t const old_size)
    {
        // if constexpr (ParticleArray::layout_mode == core::LayoutMode::AoSMapped)
        //     destParticles.map_particles(old_size);

        core::ParticleArrayService::sync<2, ParticleType_v>(destParticles);
    }

    void operator()(SAMRAI::pdat::CellOverlap const& destFieldOverlap)
    {
        // We get the source box that contains ghost region in order to get local index later
        // same for destinationGhostBox and destinationDomainBox the later will allow to get an
        // index relative to the interior
        forBoxes(destFieldOverlap.getDestinationBoxContainer());
    }

    auto static getSplitBox(core::Box<int, dim> destinationBox)
    {
        return destinationBox.grow(Splitter::maxCellDistanceFromSplit());
    }

    ~ParticlesRefining() {}
};


} // namespace PHARE::amr

namespace PHARE
{
namespace amr
{



    template<typename ParticleArray, ParticlesDataSplitType splitType, typename Splitter>
    class ParticlesRefineOperator : public SAMRAI::hier::RefineOperator
    {
    public:
        static constexpr auto dim           = Splitter::dimension;
        static constexpr auto interpOrder   = Splitter::interp_order;
        static constexpr auto nbRefinedPart = Splitter::nbRefinedPart;

        ParticlesRefineOperator()
            : SAMRAI::hier::RefineOperator{"ParticlesDataSplit_" + splitName_(splitType)}
        {
        }

        virtual ~ParticlesRefineOperator() = default;

        /** @brief a priority of 0 means that this operator
         * will be applied first
         */
        virtual int getOperatorPriority() const override { return 0; }

        virtual SAMRAI::hier::IntVector
        getStencilWidth(SAMRAI::tbox::Dimension const& dimension) const override
        {
            return SAMRAI::hier::IntVector{dimension, ghostWidthForParticles<interpOrder>()};
        }

        /** @brief perform a split and keep those that are inside a fineOverlap
         *
         */
        virtual void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                            int const destinationComponent, int const sourceComponent,
                            SAMRAI::hier::BoxOverlap const& fineOverlap,
                            SAMRAI::hier::IntVector const& /*ratio*/) const override
        {
            // For the particles we index them as a CellIndex (for the iCell)
            // ie the particles in the iCell live between lower left node of the iCell
            // and upper right of the same iCell
            auto const& destinationFieldOverlap
                = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(fineOverlap);

            // We then need to get our ParticlesData from the patch
            auto destinationParticlesData = std::dynamic_pointer_cast<ParticlesData<ParticleArray>>(
                destination.getPatchData(destinationComponent));

            auto const sourceParticlesData
                = std::dynamic_pointer_cast<ParticlesData<ParticleArray>>(
                    source.getPatchData(sourceComponent));

            ParticlesRefining<ParticleArray, splitType, Splitter>{
                *sourceParticlesData, *destinationParticlesData}(destinationFieldOverlap);
        }

    private:
        std::string splitName_(ParticlesDataSplitType splitTypeUsed)
        {
            switch (splitTypeUsed)
            {
                case ParticlesDataSplitType::coarseBoundary: return "coarseBoundary";
                case ParticlesDataSplitType::interior: return "interior";
                case ParticlesDataSplitType::coarseBoundaryOld: return "coarseBoundaryOld";
                case ParticlesDataSplitType::coarseBoundaryNew: return "coarseBoundaryNew";
                default: throw std::runtime_error("End of enum class possible range");
            }
        }
    };

} // namespace amr
} // namespace PHARE


namespace PHARE::amr
{
template<typename ParticleArray, typename Splitter>
struct RefinementParams
{
    using InteriorParticleRefineOp
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::interior, Splitter>;

    using CoarseToFineRefineOpOld
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::coarseBoundaryOld,
                                  Splitter>;

    using CoarseToFineRefineOpNew
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::coarseBoundaryNew,
                                  Splitter>;
};

} // namespace PHARE::amr


#endif
