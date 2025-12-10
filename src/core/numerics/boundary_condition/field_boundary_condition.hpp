#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>
#include <cstdint>

namespace PHARE::core
{
template<IsField FieldT, IsGridLayout GridLayoutT>
class IFieldBoundaryCondition
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using PhysicalQuantityT                = FieldT::physical_quantity_type;

    static_assert(std::is_same_v<PhysicalQuantityT, typename GridLayoutT::Quantity::Scalar>);

    virtual void apply(FieldT& field, Box<std::uint32_t, dimension> const& localGhostBox,
                       GridLayoutT const& gridLayout, double const& time)
        = 0;

    virtual BdryLoc::Type getLocation() const             = 0;
    virtual Direction getDirection() const                = 0;
    virtual Side getSide() const                          = 0;
    virtual PhysicalQuantityT getPhysicalQuantity() const = 0;
    virtual QtyCentering getQtyCentering() const          = 0;

    virtual ~IFieldBoundaryCondition() = default;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
