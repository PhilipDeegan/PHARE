#ifndef PHARE_LINEAR_WEIGHTER_H
#define PHARE_LINEAR_WEIGHTER_H


#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/field/field.h"
#include "core/utilities/constants.h"
#include "core/utilities/point/point.h"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <vector>


namespace PHARE
{
namespace amr
{
    /** @brief  This class calculates the distances of each fine index within a coarse cell
     * from the left-most coarse index of the same kind (dual/primal)
     */
    template<typename Float>
    class LinearWeighter
    {
    public:
        using FineIndexWeight  = std::array<Float, 2>;
        using FineIndexWeights = std::vector<FineIndexWeight>;


        LinearWeighter(core::QtyCentering centering, std::size_t ratio);

        std::vector<Float> const& getUniformDistances() const { return distFromLeftNode_; }
        FineIndexWeights const& weights() const { return weights_; }

    private:
        std::vector<Float> distFromLeftNode_;
        FineIndexWeights weights_;
    };


    template<typename Float, typename T, std::size_t... Is>
    std::array<LinearWeighter<Float>, sizeof...(Is)>
    make_weighters(const std::array<T, sizeof...(Is)>& values, SAMRAI::hier::IntVector ratio,
                   std::index_sequence<Is...>)
    {
        return {{(LinearWeighter<Float>{values[Is], static_cast<std::size_t>(ratio[Is])})...}};
    }
} // namespace amr
} // namespace PHARE

#include <algorithm>

namespace PHARE::amr
{
template<typename Float>
LinearWeighter<Float>::LinearWeighter(core::QtyCentering centering, std::size_t ratio)

{
    auto nbrPoints = ratio;
    assert(nbrPoints > 1);
    distFromLeftNode_.resize(nbrPoints);
    bool isEvenRatio   = ratio % 2 == 0;
    auto smallCellSize = 1. / ratio;

    std::iota(std::begin(distFromLeftNode_), std::end(distFromLeftNode_), 0);

    // when we are primal we have the coarse centering
    // that lie on top of a fine centered data
    if (centering == core::QtyCentering::primal)
    {
        std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                       std::begin(distFromLeftNode_),
                       [ratio](auto const& v) { return static_cast<Float>(v) / ratio; });
    }
    else
    {
        if (isEvenRatio)
        {
            auto middle = std::begin(distFromLeftNode_) + distFromLeftNode_.size() / 2;
            std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                           std::begin(distFromLeftNode_), [smallCellSize](auto const& v) {
                               return (0.5 + static_cast<Float>(v)) * smallCellSize;
                           });
            std::rotate(std::begin(distFromLeftNode_), middle, std::end(distFromLeftNode_));
        }

        else
        {
            auto middle = std::begin(distFromLeftNode_) + distFromLeftNode_.size() / 2 + 1;
            std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                           std::begin(distFromLeftNode_), [smallCellSize](auto const& v) {
                               return static_cast<Float>(v) * smallCellSize;
                           });

            std::rotate(std::begin(distFromLeftNode_), middle, std::end(distFromLeftNode_));
        }
    }


    std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                   std::back_inserter(weights_), [](auto const& d) {
                       return std::array<Float, 2>{{Float{1} - d, d}};
                   });
}

} // namespace PHARE::amr


#endif
