#ifndef DEFAULT_HYBRID_TAGGER_STRATEGY_H
#define DEFAULT_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include <cstddef>


namespace PHARE::amr
{
template<typename HybridModel>
class DefaultHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;
};

template<typename HybridModel>
void DefaultHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                   gridlayout_type const& layout, int* tags) const
{
    auto& Bx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    auto& Bz = model.state.electromag.B.getComponent(PHARE::core::Component::Z);

    auto& N = model.state.ions.density();

    double threshold = 0.1;

    // we loop on cell indexes for all qties regardless of their centering
    auto const& [start_x, start_y]
        = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    auto endCell_x = layout.nbrCells()[0] - 1;

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    auto tagsv
        = core::NdArrayView<dimension, int, int*, /*c_ordering=*/false>(tags, layout.nbrCells());

    if constexpr (dimension == 1)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::nbrGhosts() > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iCell = 0u, ix = start_x; iCell < endCell_x + oneOrZero; ++ix, ++iCell)
        {
            auto crit_by_x = (By(ix + 2) - By(ix)) / (1 + By(ix + 1) - By(ix));
            auto crit_bz_x = (Bz(ix + 2) - Bz(ix)) / (1 + Bz(ix + 1) - Bz(ix));
            auto criter    = std::max(crit_by_x, crit_bz_x);

            if (criter > threshold)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 2)
    {
        auto const& [start_y, end_y]
            = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);

        auto endCell_y = layout.nbrCells()[1] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x < endCell_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y < endCell_y; ++iy, ++iTag_y)
            {
                auto field_diff = [&](auto const& F) //
                {
                    return std::make_tuple(
                        std::abs((F(ix + 2, iy) - F(ix, iy)) / (1 + F(ix + 1, iy) - F(ix, iy))),
                        std::abs((F(ix, iy + 2) - F(ix, iy)) / (F(ix, iy + 1) - F(ix, iy) + 1)));
                };

                auto const& [Bx_x, Bx_y] = field_diff(Bx);
                auto const& [By_x, By_y] = field_diff(By);
                auto const& [Bz_x, Bz_y] = field_diff(Bz);
                auto crit                = std::max({Bx_x, Bx_y, By_x, By_y, Bz_x, Bz_y});

                if (crit > threshold)
                {
                    tagsv(iTag_x, iTag_y) = 1;
                }
                else
                {
                    tagsv(iTag_x, iTag_y) = 0;
                }
            }
        }
        for (auto iTag_x = 0u; iTag_x <= endCell_x; ++iTag_x)
        {
            tagsv(iTag_x, endCell_y) = tagsv(iTag_x, endCell_y - 1);
        }
        for (auto iTag_y = 0u; iTag_y <= endCell_y; ++iTag_y)
        {
            tagsv(endCell_x, iTag_y) = tagsv(endCell_x - 1, iTag_y);
        }
    }
}
} // namespace PHARE::amr

#endif // DEFAULT_HYBRID_TAGGER_STRATEGY_H
