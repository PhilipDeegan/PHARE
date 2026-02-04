#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>
#include <cstdint>

namespace PHARE::core
{
/**
 * @brief Interface for applying boundary conditions to (tensor) fields.
 *        This class provides a common interface for both Fields and TensorFields.
 * @tparam ScalarOrTensorFieldT The type of the field (must satisfy IsField or IsTensorField).
 * @tparam GridLayoutT The layout type of the grid (must satisfy IsGridLayout).
 */
template<typename ScalarOrTensorFieldT, IsGridLayout GridLayoutT>
    requires(IsField<ScalarOrTensorFieldT> || IsTensorField<ScalarOrTensorFieldT>)
class IFieldBoundaryCondition
{
public:
    /// Boolean flag indicating if the field is a scalar.
    static constexpr bool is_scalar   = IsField<ScalarOrTensorFieldT>;
    static constexpr size_t dimension = GridLayoutT::dimension;
    static constexpr size_t N = NumberOfComponentsSelector<ScalarOrTensorFieldT, is_scalar>::value;

    using physical_quantity_type
        = PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using field_type = FieldTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;

    virtual void apply(ScalarOrTensorFieldT& field,
                       Box<std::uint32_t, dimension> const& localGhostBox,
                       GridLayoutT const& gridLayout, double const& time)
        = 0;

    virtual BdryLoc::Type getLocation() const                    = 0;
    virtual Direction getDirection() const                       = 0;
    virtual Side getSide() const                                 = 0;
    virtual physical_quantity_type getPhysicalQuantity() const   = 0;
    virtual std::array<QtyCentering, N> getQtyCenterings() const = 0;

    virtual ~IFieldBoundaryCondition() = default;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
