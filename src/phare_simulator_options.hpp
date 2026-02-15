#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/def/phare_config.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE
{


struct SimOpts
{
    std::size_t dimension    = 1;
    std::size_t interp_order = 1;

    core::LayoutMode layout_mode = core::LayoutMode::AoSMapped;
    AllocatorMode alloc_mode     = AllocatorMode::CPU;

    std::size_t nbRefinedPart = core::defaultNbrRefinedParts(dimension, interp_order);

    auto static constexpr make(std::size_t const dim, std::size_t const interp,
                               std::size_t const nbRefinedPart)
    {
        return SimOpts{.dimension = dim, .interp_order = interp, .nbRefinedPart = nbRefinedPart};
    }
};



} // namespace PHARE

#endif // PHARE_SIMULATOR_OPTIONS_HPP
