#ifndef DEFAULT_HYBRID_TAGGER_STRATEGY_H
#define DEFAULT_HYBRID_TAGGER_STRATEGY_H

#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/types.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "hybrid_tagger_strategy.hpp"
#include "initializer/data_provider.hpp"

#include <cstddef>

namespace PHARE::amr
{

template<typename HybridModel>
struct DefaultHybridFieldTagger
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

    void tag(auto const& layout, int* tags, auto const& B, auto const& layout1,
             auto const& offset) const;

    double const threshold_ = 0.1;
};


template<typename HybridModel>
class DefaultHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using ParticleArray_t           = typename HybridModel::particle_array_type;
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;


public:
    DefaultHybridTaggerStrategy(initializer::PHAREDict const& dict)
        : threshold_{cppdict::get_value(dict, "threshold", 0.1)}
    {
    }
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;

private:
    double threshold_ = 0.1;
};


template<typename HybridModel>
void DefaultHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                   gridlayout_type const& layout, int* tags) const
{
    DefaultHybridFieldTagger<HybridModel> tagger{threshold_};

    auto& B = model.state.electromag.B;

    using enum core::LayoutMode;

    if constexpr (core::any_in(ParticleArray_t::layout_mode, AoSTS))
    {
        using Field_vt       = HybridModel::field_type::value_type;
        using TensorField_vt = core::basic::TensorField<Field_vt, 1>;

        auto const ntiles = B[0]().size();
        for (std::size_t tidx = 0; tidx < ntiles; ++tidx)
        {
            auto Btile = B.template as<TensorField_vt>([&](auto& c) { return c()[tidx]; });

            auto const& tile_layout   = B[0]()[tidx].layout();
            auto const tag_idx_offset = layout.AMRToLocal(B[0]()[tidx].layout().AMRBox().lower);

            tagger.tag(layout, tags, Btile, tile_layout,
                       tag_idx_offset - gridlayout_type::nbrGhosts());
        }
    }
    else
    {
        tagger.tag(layout, tags, B, layout, core::Point{ConstArray<int, dimension>()});
    }
}

template<typename HybridModel>
void DefaultHybridFieldTagger<HybridModel>::tag(auto const& layout, int* tags, auto const& B,
                                                auto const& layout1, auto const& offset) const
{
    auto& [Bx, By, Bz] = B();

    // we loop on cell indexes for all qties regardless of their centering
    // auto const& [start[0], _]
    //     = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    auto const domain_box = layout1.domainBoxFor(core::QtyCentering::dual);
    auto const start      = domain_box.lower;
    auto const end        = core::Point{for_N_make_array<dimension>(
        [&](auto i) { return layout1.nbrCells()[i] - 1; })}; // domain_box.upper - 1;

    // override end[0] because the tag buffer does not have ghost cells
    // and physicalEnd will account for ghost cells
    // auto const& end[0] = layout.nbrCells()[0] - 1;

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    bool constexpr c_ordering = false;
    auto tagsv = core::NdArrayView<dimension, int, c_ordering>(tags, layout.nbrCells());

    if constexpr (dimension == 1 and false)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        for (auto iCell = 0u, ix = start[0]; iCell <= end[0]; ++ix, ++iCell)
        {
            auto crit_by_x = (By(ix + 2) - By(ix)) / (1 + By(ix + 1) - By(ix));
            auto crit_bz_x = (Bz(ix + 2) - Bz(ix)) / (1 + Bz(ix + 1) - Bz(ix));
            auto criter    = std::max(crit_by_x, crit_bz_x);

            if (criter > threshold_)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 1)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::nbrGhosts() > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;

        for (auto iCell = 0u, ix = start[0]; iCell < end[0] + oneOrZero; ++ix, ++iCell)
        {
            auto Byavg     = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2));
            auto Bzavg     = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2));
            auto Byavgp1   = 0.2 * (By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3));
            auto Bzavgp1   = 0.2 * (Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3));
            auto criter_by = std::abs(Byavgp1 - Byavg) / (1 + std::abs(Byavg));
            auto criter_bz = std::abs(Bzavgp1 - Bzavg) / (1 + std::abs(Bzavg));
            auto criter_b  = std::sqrt(criter_by * criter_by + criter_bz * criter_bz);
            auto criter    = criter_b;

            auto const tag_idx = core::Point{iCell} + offset;
            if (criter > threshold_)
            {
                tagsv(tag_idx) = 1;
            }
            else
                tagsv(tag_idx) = 0;
        }
    }
    if constexpr (dimension == 2)
    {
        for (auto iTag_x = 0u, ix = start[0]; iTag_x <= end[0]; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start[1]; iTag_y <= end[1]; ++iy, ++iTag_y)
            {
                auto field_diff = [&](auto const& F) //
                {
                    auto const delta_2x = std::abs(F(ix + 2, iy) - F(ix, iy));
                    auto const delta_2y = std::abs(F(ix, iy + 2) - F(ix, iy));
                    auto const delta_x  = std::abs(F(ix + 1, iy) - F(ix, iy));
                    auto const delta_y  = std::abs(F(ix, iy + 1) - F(ix, iy));

                    auto const criter_x = delta_2x / (1 + delta_x);
                    auto const criter_y = delta_2y / (1 + delta_y);

                    return std::make_tuple(criter_x, criter_y);
                };

                auto const& [Bx_x, Bx_y] = field_diff(Bx);
                auto const& [By_x, By_y] = field_diff(By);
                auto const& [Bz_x, Bz_y] = field_diff(Bz);
                auto crit                = std::max({Bx_x, Bx_y, By_x, By_y, Bz_x, Bz_y});

                auto const tag_idx = core::Point{iTag_x, iTag_y} + offset;
                if (crit > threshold_)
                {
                    tagsv(tag_idx) = 1;
                }
                else
                {
                    tagsv(tag_idx) = 0;
                }
            }
        }
    }
}

} // namespace PHARE::amr

#endif // DEFAULT_HYBRID_TAGGER_STRATEGY_H
