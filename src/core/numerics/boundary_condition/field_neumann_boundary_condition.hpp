#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP

#include "core/numerics/boundary_condition/field_boundary_condition_dispatcher.hpp"

namespace PHARE::core
{
template<IsField FieldT, IsGridLayout GridLayoutT>
class FieldNeumannBoundaryCondition
    : public FieldBoundaryConditionDispatcher<FieldT, GridLayoutT,
                                              FieldNeumannBoundaryCondition<FieldT, GridLayoutT>>
{
public:
    using PhysicalQuantityT = FieldT::physical_quantity_type;
    using Super
        = FieldBoundaryConditionDispatcher<FieldT, GridLayoutT,
                                           FieldNeumannBoundaryCondition<FieldT, GridLayoutT>>;

    using Super::dimension;

    FieldNeumannBoundaryCondition() = delete;
    FieldNeumannBoundaryCondition(BdryLoc::Type const& location,
                                  PhysicalQuantityT const& physicalQuantity)
        : Super(location, physicalQuantity)
    {
    } //

    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition const&)            = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition const&) = default;
    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition&&)                 = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition&&)      = default;

    virtual ~FieldNeumannBoundaryCondition() = default;

    template<Direction direction, Side side, QtyCentering centering>
    void apply_specialized(FieldT& field, Box<std::uint32_t, dimension> const& localGhostBox,
                           GridLayoutT const& gridLayout, double const time)
    {
        using Index = Point<std::uint32_t, dimension>;

        Index physicalLimitIndex = (side == Side::LOWER) ? gridLayout.physicalStartIndex(centering)
                                                         : gridLayout.physicalEndIndex(centering);
        for (Index const& index : localGhostBox)
        {
            Index mirrorIndex = boundary_mirrored<dimension, direction, side, centering>(
                index, physicalLimitIndex);
            field(index) = field(mirrorIndex);
        };
    }
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
