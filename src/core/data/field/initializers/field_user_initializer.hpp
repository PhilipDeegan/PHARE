#ifndef _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_

#include "core/vector.hpp"

#include <tuple>

namespace PHARE::core
{

class FieldUserFunctionInitializer
{
public:
    template<typename Field, typename GridLayout, typename InitFunction>
    void static initialize(Field& field, GridLayout const& layout, InitFunction const& init)
    {
        auto constexpr static alloc_mode = Field::alloc_mode;
        // static_assert(alloc_mode == AllocatorMode::GPU_UNIFIED);

        auto const indices = layout.ghostStartToEndIndices(field, /*includeEnd=*/true);
        auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
            indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
            });

        // keep grid data alive
        auto grid = std::apply([&](auto const&... args) { return init(args...); }, coords);
        assert(field.size() == grid->size());

        if constexpr (alloc_mode == AllocatorMode::CPU || alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
                std::apply([&](auto&... args) { field(args...) = (*grid)[cell_idx]; },
                           indices[cell_idx]);
        }
        else
        {
            if constexpr (CompileOptions::WithUmpire and CompileOptions::WithRAJA)
            {
                PHARE_WITH_RAJA(PHARE::core::raja::copy(field.data(), grid->data(), field.size()));
            }
            else
            {
                throw std::runtime_error("PHARE::core::initialize_field NO ALTERNATIVE");
            }
        }
    }
};


} // namespace PHARE::core

#endif
