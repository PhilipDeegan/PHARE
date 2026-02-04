#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_dispatcher.hpp"

#include <cstddef>
#include <tuple>

namespace PHARE::core
{
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldNeumannBoundaryCondition
    : public FieldBoundaryConditionDispatcher<
          ScalarOrTensorFieldT, GridLayoutT,
          FieldNeumannBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>
{
public:
    using Super = FieldBoundaryConditionDispatcher<
        ScalarOrTensorFieldT, GridLayoutT,
        FieldNeumannBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>;
    using physical_quantity_type = Super::physical_quantity_type;
    using field_type             = Super::field_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr size_t N         = Super::N;
    static constexpr bool is_scalar   = Super::is_scalar;

    FieldNeumannBoundaryCondition() = delete;
    FieldNeumannBoundaryCondition(BdryLoc::Type const& location,
                                  physical_quantity_type const& physicalQuantity)
        : Super(location, physicalQuantity)
    {
    } //

    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition const&)            = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition const&) = default;
    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition&&)                 = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition&&)      = default;

    virtual ~FieldNeumannBoundaryCondition() = default;

    template<Direction direction, Side side, QtyCentering... Centerings>
    void apply_specialized(ScalarOrTensorFieldT& scalarOrTensorField,
                           Box<std::uint32_t, dimension> const& localGhostBox,
                           GridLayoutT const& gridLayout, double const time)
    {
        using Index = Point<std::uint32_t, dimension>;

        constexpr std::array<QtyCentering, N> centerings = {Centerings...};

        // no other way than using a lambda builder
        auto fields = [&]() {
            if constexpr (is_scalar)
                return std::make_tuple(scalarOrTensorField);
            else
                return scalarOrTensorField.components();
        }();

        for_N<N>([&](auto i) {
            constexpr QtyCentering centering = centerings[i];
            field_type& field                = std::get<i>(fields);
            auto fieldBox = gridLayout.toFieldBox(localGhostBox, field.physicalQuantity());
            Index physicalLimitIndex = (side == Side::LOWER)
                                           ? gridLayout.physicalStartIndex(centering)
                                           : gridLayout.physicalEndIndex(centering);
            for (Index const& index : fieldBox)
            {
                Index mirrorIndex = boundary_mirrored<dimension, direction, side, centering>(
                    index, physicalLimitIndex);
                field(index) = field(mirrorIndex);
            }
        });
    }
}; // class FieldNeumannBoundaryCondition

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
