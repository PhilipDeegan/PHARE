#ifndef PHARE_LINEAR_WEIGHTER_HPP
#define PHARE_LINEAR_WEIGHTER_HPP


#include "core/def/phare_mpi.hpp"


#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"

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
    class LinearWeighter
    {
    public:
        using FineIndexWeight  = std::array<core::floater_t<0>, 2>;
        using FineIndexWeights = std::vector<FineIndexWeight>;


        LinearWeighter(core::QtyCentering centering, std::size_t ratio);

        std::vector<core::floater_t<4>> const& getUniformDistances() const
        {
            return distFromLeftNode_;
        }
        FineIndexWeights const& weights() const { return weights_; }

    private:
        std::vector<core::floater_t<4>> distFromLeftNode_;
        FineIndexWeights weights_;
    };


    template<typename T, std::size_t... Is>
    NO_DISCARD std::array<LinearWeighter, sizeof...(Is)>
    make_weighters(const std::array<T, sizeof...(Is)>& values, SAMRAI::hier::IntVector ratio,
                   std::index_sequence<Is...>)
    {
        return {{(LinearWeighter{values[Is], static_cast<std::size_t>(ratio[Is])})...}};
    }
} // namespace amr
} // namespace PHARE


#endif
