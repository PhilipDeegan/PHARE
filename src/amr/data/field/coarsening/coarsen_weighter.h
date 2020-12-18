#ifndef PHARE_COARSEN_WEIGHTER_H
#define PHARE_COARSEN_WEIGHTER_H

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>


namespace PHARE
{
namespace amr
{
    /**
     * @brief The CoarsenWeighter class computes the weights to use for each nbrPoints fine nodes to
     * get the value on the associated coarse node
     */
    template<typename Float>
    class CoarsenWeighter
    {
    public:
        explicit CoarsenWeighter(std::size_t nbrPoints)
        {
            assert(nbrPoints > 1); // we want to have at least 2 points for coarsening operations
            computeWeights_(nbrPoints);
        }

        std::vector<Float> const& weights() const { return weights_; }

    private:
        std::vector<Float> weights_;

        Float findX_(std::size_t nbrPoints) const;
        void computeWeights_(std::size_t nbrPoints);
    };



    template<typename Float, typename T, std::size_t... Is>
    std::array<CoarsenWeighter<Float>, sizeof...(Is)>
    make_weighters(const std::array<T, sizeof...(Is)>& nbrPoints, std::index_sequence<Is...>)
    {
        return {{(CoarsenWeighter<Float>{nbrPoints[Is]})...}};
    }
} // namespace amr

} // namespace PHARE



namespace PHARE::amr
{
template<typename Float>
Float CoarsenWeighter<Float>::findX_(std::size_t nbrPoints) const
{
    Float x = 0.;

    if (nbrPoints % 2 != 0)
    {
        x = 1.;
        for (std::size_t i = 1; i <= (nbrPoints - 1) / 2; ++i)
        {
            x += 2 * 1. / (i + 1);
        }
    }
    else
    {
        for (std::size_t i = 1; i <= nbrPoints / 2; ++i)
        {
            x += 2 * 1. / i;
        }
    }

    return x;
}




template<typename Float>
void CoarsenWeighter<Float>::computeWeights_(std::size_t nbrPoints)
{
    weights_.resize(nbrPoints);

    auto x = findX_(nbrPoints);


    if (nbrPoints % 2 != 0)
    {
        auto const halfIndex = (nbrPoints - 1) / 2;

        auto const halfNumberOfPointsLeft
            = nbrPoints / 2; // half of the points needed besides the one on the middle

        weights_[halfIndex] = 1. / x;

        for (std::size_t i = 1; i <= halfNumberOfPointsLeft; ++i)
        {
            Float factor            = static_cast<Float>(i + 1);
            weights_[halfIndex - i] = 1. / (factor * x);
            weights_[halfIndex + i] = 1. / (factor * x);
        }
    }




    else
    {
        auto const halfIndexRight = nbrPoints / 2;
        auto const halfIndexLeft  = halfIndexRight - 1;

        weights_[halfIndexRight] = 1. / x;
        weights_[halfIndexLeft]  = 1. / x;

        auto const halfNumberOfPointsLeft = (nbrPoints / 2) - 1;

        for (std::size_t i = 1; i <= halfNumberOfPointsLeft; ++i)
        {
            Float factor = static_cast<Float>(i + 1);

            weights_[halfIndexRight + i] = 1. / (factor * x);
            weights_[halfIndexLeft - i]  = 1. / (factor * x);
        }
    }
}

} // namespace PHARE::amr


#endif
