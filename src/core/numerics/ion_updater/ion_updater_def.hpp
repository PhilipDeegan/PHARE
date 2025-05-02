#ifndef PHARE_ION_UPDATER_DEF_HPP
#define PHARE_ION_UPDATER_DEF_HPP


#include "core/data/particles/particle_array_def.hpp"
#include "core/utilities/box/box.hpp"

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
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;

    GridLayout_t const layout;
    Box_t const nonLevelGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);
};


template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;

    UpdaterCellMapSelectionBoxing(auto&&... args)
        : Super{args...}
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
};


} // namespace PHARE::core


// ###############################################

#if PHARE_HAVE_MKN_GPU || defined(_MKN_WITH_MKN_KUL_)
#include "mkn/kul/os.hpp"

namespace PHARE::core::detail
{
auto static const timings_dir_str
    = get_env_as("PHARE_ASYNC_TIMES", std::string{".phare/async/multi_updater"});

static bool ion_updater_io_setup = []() {
    PHARE_LOG_LINE_SS("MAKE DIR!")
    mkn::kul::Dir timings{timings_dir_str};
    timings.mk();
    return true;
}();

} // namespace PHARE::core::detail

#endif // PHARE_HAVE_MKN_GPU

// ###############################################


#endif // ION_UPDATER_DEF_HPP
