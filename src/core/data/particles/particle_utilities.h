#ifndef PHARE_PARTICLE_UTILITIES
#define PHARE_PARTICLE_UTILITIES

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/particles/particle.h"
#include "core/utilities/point/point.h"

#include <array>


namespace PHARE
{
namespace core
{
    /**
     * @brief positionAsPoint returns a point holding the physical position of the macroparticle.
     * The function assumes the iCell of the particle is in AMR index space.
     */
    template<typename Particle, typename GridLayout>
    decltype(auto) positionAsPoint(Particle const& particle, GridLayout const& layout)
    {
        static_assert(Particle::dimension == GridLayout::dimension,
                      "Invalid use of function, dimension template mismatch");

        using Float = typename Particle::float_type;

        Point<Float, GridLayout::dimension> position;
        auto origin       = layout.origin();
        auto startIndexes = layout.physicalStartIndex(QtyCentering::primal);
        auto meshSize     = layout.meshSize();
        auto iCell        = layout.AMRToLocal(Point{particle.iCell});

        for (auto iDim = 0u; iDim < GridLayout::dimension; ++iDim)
        {
            auto delta     = static_cast<Float>(particle.delta[iDim]);
            position[iDim] = origin[iDim];
            position[iDim] += (iCell[iDim] - startIndexes[iDim] + delta) * meshSize[iDim];
        }
        return position;
    }
} // namespace core




} // namespace PHARE



#endif
