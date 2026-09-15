#ifndef PHARE_ION_UPDATER_DEF_HPP
#define PHARE_ION_UPDATER_DEF_HPP


#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstdint>

namespace PHARE::core
{
enum class UpdaterMode : std::uint16_t { domain_only = 0, all };

template<typename GridLayout>
struct UpdaterSelectionBoxing;

template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing;

template<typename Updater_t, typename GridLayout>
constexpr auto* cellmapselection_boxing_impl()
{
    using Selector_t = Updater_t::Pusher::ParticleSelector;
    return static_cast<UpdaterCellMapSelectionBoxing<Selector_t, GridLayout>*>(0);
}

template<typename Updater_t, typename GridLayout>
constexpr auto* selection_boxing_impl()
{
    using enum LayoutMode;
    if constexpr (Updater_t::ParticleArray_t::layout_mode == AoSMapped)
        return cellmapselection_boxing_impl<Updater_t, GridLayout>();
    else
        return static_cast<UpdaterSelectionBoxing<GridLayout>*>(0);
}



template<typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;

    UpdaterSelectionBoxing(GridLayout_t const& layout_, std::vector<Box_t> const& nonLevelGhostBox_)
        : layout{layout_}
        , nonLevelGhostBox{nonLevelGhostBox_}
    {
    }

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);
};


template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;

    UpdaterCellMapSelectionBoxing(GridLayout_t const& layout_,
                                  std::vector<Box_t> const& nonLevelGhostBox_)
        : Super{layout_, nonLevelGhostBox_}
    {
    }

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inDomainBox = [domainBox = Super::domainBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return core::isIn(cell, domainBox); });
    };

    Selector_t const inGhostBox = [ghostBox = Super::ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return isIn(cell, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [nonLevelGhostBox = Super::nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(
                  particleRange, [&](auto const& cell) { return isIn(cell, nonLevelGhostBox); });
          };

    Selector_t const inGhostLayer
        = [ghostBox = Super::ghostBox, domainBox = Super::domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, domainBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = Super::ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return !isIn(cell, ghostBox); });
    };
};



template<typename GridLayout>
struct UpdaterTileSetSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static dimension      = GridLayout::dimension;
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;

    UpdaterTileSetSelectionBoxing(GridLayout_t const& layout_,
                                  std::vector<Box_t> const& nonLevelGhostBox_)
        : Super{layout_, nonLevelGhostBox_}
    {
    }
};


} // namespace PHARE::core


#endif // ION_UPDATER_DEF_HPP
