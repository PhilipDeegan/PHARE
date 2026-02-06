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
 * @brief Interface for applying boundary conditions to scalar or tensor fields.
 *
 * Concrete field boundary conditions are provided by implementating this interface.
 *
 * @tparam ScalarOrTensorFieldT The type of the scalarOrTensorField (must satisfy IsField or
 * IsTensorField).
 * @tparam The grid layout type (must satisfy IsGridLayout).
 *
 */
template<typename ScalarOrTensorFieldT, IsGridLayout GridLayoutT>
    requires(IsField<ScalarOrTensorFieldT> || IsTensorField<ScalarOrTensorFieldT>)
class IFieldBoundaryCondition
{
public:
    static constexpr bool is_scalar   = IsField<ScalarOrTensorFieldT>;
    static constexpr size_t dimension = GridLayoutT::dimension;
    static constexpr size_t N = NumberOfComponentsSelector<ScalarOrTensorFieldT, is_scalar>::value;

    using physical_quantity_type
        = PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using field_type = FieldTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;

    /**
     * @brief Enforce the boundary condition on the provided scalar/tensor @p scalarOrTensorField,
     * by filling accordingly the ghost cells contained in the local box @p localGhostBox, at the
     * physical time @p time, and considering that the boundary is located at @p boundaryLocation.
     *
     * @param scalarOrTensorField The scalar or tensor to which we apply the boundary condition.
     * @param boundaryLocation The location of the physical boundary.
     * @param localGhostBox The box containing the ghost cells/nodes to fill.
     * @param gridLayout The grid layout.
     * @param time The physical time, useful for time-dependant boundary conditions.
     *
     */
    virtual void apply(ScalarOrTensorFieldT& scalarOrTensorField,
                       BdryLoc::Type const boundaryLocation,
                       Box<std::uint32_t, dimension> const& localGhostBox,
                       GridLayoutT const& gridLayout, double const time)
        = 0;

    virtual ~IFieldBoundaryCondition() = default;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
