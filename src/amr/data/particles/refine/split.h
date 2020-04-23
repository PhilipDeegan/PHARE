#ifndef PHARE_SPLIT_H
#define PHARE_SPLIT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "core/utilities/types.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/particles/particle.h"
#include "core/utilities/box/box.h"

#include "kul/log.hpp"

namespace PHARE::amr
{
template<size_t _dimension, size_t _refinedParticlesNbr>
class ASplitter
{
public:
    static constexpr size_t dimension           = _dimension;
    static constexpr size_t refinedParticlesNbr = _refinedParticlesNbr;
    static constexpr size_t refinementFactor    = 2;

    inline void operator()(core::Particle<dimension> const& coarsePartOnRefinedGrid,
                           std::vector<core::Particle<dimension>>& refinedParticles) const
    {
        auto get = [&](size_t rpIndex, size_t index, auto& icell, auto& delta) {
            delta[index] += deltas_[rpIndex][index] * refinementFactor;
            float integra = std::floor(delta[index]);
            delta[index] -= integra;
            icell[index] += static_cast<int32_t>(integra);
        };

        for (uint32_t refinedParticleIndex = 0; refinedParticleIndex < refinedParticlesNbr;
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

protected:
    ASplitter()
        : weights_(refinedParticlesNbr)
        , deltas_(refinedParticlesNbr)
    {
    }

    std::vector<float> weights_;
    std::vector<std::array<uint32_t, dimension>> iCells_;
    std::vector<std::array<float, dimension>> deltas_;
};
} // namespace PHARE::amr

#include "split_1d.h"
#include "split_2d.h"
#include "split_3d.h"

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


    template<std::size_t dim, std::size_t interp, std::size_t nbRefinedPart>
    using SuperSplitter = typename std::conditional<
        dim == 1, Splitter_1d<interp, _nbRefinedPart>,
        typename std::conditional<dim == 2, Splitter_2d<interp, _nbRefinedPart>,
                                  Splitter_3d<interp, _nbRefinedPart>>::type>::type;

    template<std::size_t _dimension, std::size_t _interpOrder, std::size_t _nbRefinedPart>
    class Splitter : public SuperSplitter<_dimension, _interpOrder, _nbRefinedPart>
    {
    public:
        using Super = SuperSplitter<_dimension, _interpOrder, _nbRefinedPart>;

        static constexpr size_t dimension     = _dimension;
        static constexpr size_t interpOrder   = _interpOrder;
        static constexpr size_t nbRefinedPart = _nbRefinedPart;

        // to be removed
        Splitter(core::Point<int32, dimension>, uint32)
            : Super{}
        {
        }

        static constexpr int maxCellDistanceFromSplit()
        {
            return std::ceil((interpOrder + 1) * 0.5);
        }
    };

    template<std::size_t dimension, std::size_t interpOrder>
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
        constexpr static std::array<float, 3> tabD1RF2N02Delta_  = {{0.277f, 0.332f, 0.376f}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 3
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N03Weight_
            = {{{{0.5f, 0.468f, 0.474f}}, {{0.25f, 0.266f, 0.263f}}}};
        constexpr static std::array<float, 3> tabD1RF2N03Delta_ = {{0.5f, 0.556f, 0.638f}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 4
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Weight_
            = {{{{0.0f, 0.125f, 0.135f}}, {{0.0f, 0.375f, 0.365f}}}};
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Delta_
            = {{{{0.0f, 0.75f, 0.833f}}, {{0.0f, 0.25f, 0.272f}}}};


        // dimension = 1, refinement factor = 2, nbrOfBabies = 5
        constexpr static std::array<std::array<float, 3>, 3> tabD1RF2N05Weight_
            = {{{{0.0f, 0.0f, 0.375f}}, {{0.0f, 0.0f, 0.0625f}}, {{0.0f, 0.0f, 0.25f}}}};
        constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N05Delta_
            = {{{{0.0f, 0.0f, 1.f}}, {{0.0f, 0.0f, 0.5f}}}};

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
