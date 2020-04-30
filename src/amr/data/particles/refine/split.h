#ifndef PHARE_SPLIT_H
#define PHARE_SPLIT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "core/def.h" // GPU STUFF

#include "core/utilities/types.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/particles/particle.h"
#include "core/utilities/box/box.h"

namespace PHARE::amr
{
// see meta_utilities.h for list of declared permutations.
template<size_t _dimension, size_t _interp_order, size_t _nbRefinedPart>
class ASplitter
{
public:
    static constexpr size_t dimension     = _dimension;
    static constexpr size_t interp_order  = _interp_order;
    static constexpr size_t nbRefinedPart = _nbRefinedPart;

    inline void operator()(core::Particle<dimension> const& coarsePartOnRefinedGrid,
                           std::vector<core::Particle<dimension>>& refinedParticles) const
    {
        auto get = [&](size_t rpIndex, size_t index, auto& icell, auto& delta) {
            delta[index] += deltas_[rpIndex][index];
            float integra = std::floor(delta[index]);
            delta[index] -= integra;
            icell[index] += static_cast<int32_t>(integra);
        };

        for (uint32_t refinedParticleIndex = 0; refinedParticleIndex < nbRefinedPart;
             ++refinedParticleIndex)
        {
            std::array<int32_t, dimension> icell{coarsePartOnRefinedGrid.iCell};
            std::array<float, dimension> delta{coarsePartOnRefinedGrid.delta};

            for (size_t i = 0; i < dimension; i++)
                get(refinedParticleIndex, i, icell, delta);

            float weight = coarsePartOnRefinedGrid.weight * weights_[refinedParticleIndex];
            refinedParticles.push_back(
                {weight, coarsePartOnRefinedGrid.charge, icell, delta, coarsePartOnRefinedGrid.v});
        }
    }

    static constexpr int maxCellDistanceFromSplit() { return std::ceil((interp_order + 1) * 0.5); }

protected:
    ASplitter()
        : weights_(nbRefinedPart)
        , deltas_(nbRefinedPart)
    {
    }


    template<typename Weights, typename Deltas>
    inline ASplitter(Weights, Deltas);

    std::vector<float> weights_;
    std::vector<std::array<float, dimension>> deltas_;
};

template<std::size_t dimension, std::size_t interp_order, size_t nbRefinedPart>
class Splitter : public ASplitter<dimension, interp_order, nbRefinedPart>
{
    // Unspecialized template class, never to be instantiated
    Splitter() = delete;
};

template<size_t dim, size_t nbrRefineParts>
struct SplitInnerSetter
{
    // Unspecialized template class, never to be instantiated
    SplitInnerSetter() = delete;
};

} // namespace PHARE::amr

#include "split_1d.h"
#include "split_2d.h"
#include "split_3d.h"

namespace PHARE::amr
{
template<size_t dim, size_t io, size_t nb>
template<typename Weights, typename Deltas>
ASplitter<dim, io, nb>::ASplitter(Weights weight_vals, Deltas delta_vals)
    : ASplitter(/*set vector sizes*/)
{
    using Setter = SplitInnerSetter<dimension, nbRefinedPart>;
    Setter::set_weights(weights_, weight_vals);
    Setter::set_deltas(deltas_, delta_vals);
}
} // namespace PHARE::amr

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    using core::int32;
    using core::uint32;

#define ISNOTTABULATED(dim, RF) (dim != 1 || RF != 2)



    template<std::size_t _dimension, std::size_t _interpOrder, std::size_t _nbRefinedPart>
    class _Split
    {
    public:
        static constexpr size_t dimension     = _dimension;
        static constexpr size_t interpOrder   = _interpOrder;
        static constexpr size_t nbRefinedPart = _nbRefinedPart;

    private:
        std::vector<float> weights_;
        std::vector<uint32> iCellsX_;
        std::vector<float> deltasX_;
        core::Point<int32, dimension> refinementFactor_;


        // dimension = 1, refinement factor = 2, nbrOfBabies = 2
        constexpr static std::array<float, 3> tabD1RF2N02Weight_ = {{0.5, 0.5, 0.5}};
        constexpr static std::array<float, 3> tabD1RF2N02Delta_
            = {{0.551569f, 0.663959f, 0.752399f}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 3
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N03Weight_
            = {{{{0.5f, 0.468137f, 0.473943f}}, {{0.25f, 0.265931f, 0.263028f}}}};
        constexpr static std::array<float, 3> tabD1RF2N03Delta_ = {{1.0f, 1.112033f, 1.275922f}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 4
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Weight_
            = {{{{0.0f, 0.375f, 0.364766f}}, {{0.0f, 0.125f, 0.135234f}}}};
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Delta_
            = {{{{0.0f, 0.5f, 0.542949f}}, {{0.0f, 1.5f, 1.664886f}}}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 5
        constexpr static std::array<std::array<float, 3>, 3> tabD1RF2N05Weight_
            = {{{{0.0f, 0.0f, 0.375f}}, {{0.0f, 0.0f, 0.25f}}, {{0.0f, 0.0f, 0.0625f}}}};
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N05Delta_
            = {{{{0.0f, 0.0f, 1.0f}}, {{0.0f, 0.0f, 2.0f}}}};

        const std::array<std::vector<int>, 3> tabNbrOfBabies_ = {{{2, 3}, {2, 3, 4}, {2, 3, 4, 5}}};




    public:
        _Split(core::Point<int32, dimension> refineFactor)
            : refinementFactor_{refineFactor}
        {
            deltasX_.assign(nbRefinedPart, 0);
            weights_.assign(nbRefinedPart, 0);

            std::vector<int> const& validRefinedParticleNbr = tabNbrOfBabies_[interpOrder - 1];

            for (auto iDim = 0u; iDim < dimension; ++iDim)
            {
                if (ISNOTTABULATED(dimension, refineFactor[iDim]))
                {
                    std::cout << "dimension and/or refinement factor not tabulated" << std::endl;
                }

                else
                {
                    if (std::find(std::begin(validRefinedParticleNbr),
                                  std::end(validRefinedParticleNbr), nbRefinedPart)
                        == std::end(validRefinedParticleNbr))
                    {
                        throw std::runtime_error("Invalid refined particle number");
                    }

                    else
                    {
                        // weights & deltas are coming from the tabulated values
                        switch (nbRefinedPart)
                        {
                            case 2:
                                weights_[0] = tabD1RF2N02Weight_[interpOrder - 1];
                                weights_[1] = tabD1RF2N02Weight_[interpOrder - 1];

                                deltasX_[0] = -tabD1RF2N02Delta_[interpOrder - 1];
                                deltasX_[1] = +tabD1RF2N02Delta_[interpOrder - 1];
                                break;

                            case 3:
                                weights_[0] = tabD1RF2N03Weight_[0][interpOrder - 1];
                                weights_[1] = tabD1RF2N03Weight_[1][interpOrder - 1];
                                weights_[2] = tabD1RF2N03Weight_[1][interpOrder - 1];

                                deltasX_[0] = 0.0;
                                deltasX_[1] = -tabD1RF2N03Delta_[interpOrder - 1];
                                deltasX_[2] = +tabD1RF2N03Delta_[interpOrder - 1];
                                break;

                            case 4:
                                weights_[0] = tabD1RF2N04Weight_[0][interpOrder - 1];
                                weights_[1] = tabD1RF2N04Weight_[0][interpOrder - 1];
                                weights_[2] = tabD1RF2N04Weight_[1][interpOrder - 1];
                                weights_[3] = tabD1RF2N04Weight_[1][interpOrder - 1];

                                deltasX_[0] = -tabD1RF2N04Delta_[0][interpOrder - 1];
                                deltasX_[1] = +tabD1RF2N04Delta_[0][interpOrder - 1];
                                deltasX_[2] = -tabD1RF2N04Delta_[1][interpOrder - 1];
                                deltasX_[3] = +tabD1RF2N04Delta_[1][interpOrder - 1];
                                break;

                            case 5:
                                weights_[0] = tabD1RF2N05Weight_[0][interpOrder - 1];
                                weights_[1] = tabD1RF2N05Weight_[1][interpOrder - 1];
                                weights_[2] = tabD1RF2N05Weight_[1][interpOrder - 1];
                                weights_[3] = tabD1RF2N05Weight_[2][interpOrder - 1];
                                weights_[4] = tabD1RF2N05Weight_[2][interpOrder - 1];

                                deltasX_[0] = 0.0;
                                deltasX_[1] = -tabD1RF2N05Delta_[0][interpOrder - 1];
                                deltasX_[2] = +tabD1RF2N05Delta_[0][interpOrder - 1];
                                deltasX_[3] = -tabD1RF2N05Delta_[1][interpOrder - 1];
                                deltasX_[4] = +tabD1RF2N05Delta_[1][interpOrder - 1];
                                break;

                            default: throw std::runtime_error("Invalid refined particle number");
                        }
                    }
                }
            }
        }

        ~_Split() = default;


        static constexpr int maxCellDistanceFromSplit()
        {
            return std::ceil((interpOrder + 1) * 0.5);
        }




        inline void operator()(core::Particle<dimension> const& coarsePartOnRefinedGrid,
                               std::vector<core::Particle<dimension>>& refinedParticles) const
        {
            for (uint32 refinedParticleIndex = 0; refinedParticleIndex < nbRefinedPart;
                 ++refinedParticleIndex)
            {
                if constexpr (dimension == 1)
                {
                    // the values for icell & delta are only working for 1 dim...
                    float weight = coarsePartOnRefinedGrid.weight * weights_[refinedParticleIndex];
                    int32 icell  = coarsePartOnRefinedGrid.iCell[0];
                    float delta  = coarsePartOnRefinedGrid.delta[0]
                                  + deltasX_[refinedParticleIndex] * refinementFactor_[dirX];

                    // weights & deltas are the only known values for the babies.
                    // so the icell values of each baby needs to be calculated
                    float integra = std::floor(delta);
                    delta -= integra;
                    icell += static_cast<int32>(integra);

                    refinedParticles.push_back({weight,
                                                coarsePartOnRefinedGrid.charge,
                                                {{icell}},
                                                {{delta}},
                                                coarsePartOnRefinedGrid.v});
                }
                else // unsupported dimension
                {
                    static_assert("Only 1D is supported for split at the moment");
                }
            }
        }
    };

} // namespace amr


} // namespace PHARE
#endif // endif SPLIT_H
