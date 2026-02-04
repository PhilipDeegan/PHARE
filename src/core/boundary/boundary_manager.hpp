#ifndef PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER
#define PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"

#include "initializer/data_provider.hpp"

#include <concepts>
#include <map>
#include <memory>
#include <initializer_list>
#include <stdexcept>

namespace PHARE::core
{

template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT>
class BoundaryManager
{
public:
    using scalar_quantity_type = FieldT::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename PhysicalQuantityT::Scalar>);
    using vector_field_type     = VecField<FieldT, PhysicalQuantityT>;
    using scalar_condition_type = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using vector_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;

    BoundaryManager() = delete;
    BoundaryManager(PHARE::initializer::PHAREDict const& dict,
                    std::vector<typename PhysicalQuantityT::Scalar> const& scalarQuantities,
                    std::vector<typename PhysicalQuantityT::Vector> const& vectorQuantities)
    {
        // do not actually read anything, and set everyone to Neumann for the moment
        auto locations = {BdryLoc::Type::XLO, BdryLoc::Type::XHI};
        for (auto const location : locations)
        {
            for (auto const quantity : scalarQuantities)
            {
                // we will need a factory here
                scalar_conditions_[_scalar_key_type({location, quantity})]
                    = std::make_shared<FieldNeumannBoundaryCondition<FieldT, GridLayoutT>>(
                        location, quantity);
            };

            for (auto const quantity : vectorQuantities)
            {
                // we will need a factory here
                vector_conditions_[_vector_key_type({location, quantity})] = std::make_shared<
                    FieldNeumannBoundaryCondition<vector_field_type, GridLayoutT>>(location,
                                                                                   quantity);
            };
        };
    }

    template<typename TensorPhysicalQuantity>
    auto getBoundaryCondition(BdryLoc::Type location, TensorPhysicalQuantity quantity) const
    {
        if constexpr (std::same_as<TensorPhysicalQuantity, typename PhysicalQuantityT::Scalar>)
        {
            auto it = scalar_conditions_.find({location, quantity});
            return (it != scalar_conditions_.end()) ? it->second : nullptr;
        }
        else if constexpr (std::same_as<TensorPhysicalQuantity, typename PhysicalQuantityT::Vector>)
        {
            auto it = vector_conditions_.find({location, quantity});
            return (it != vector_conditions_.end()) ? it->second : nullptr;
        }
        else
        {
            throw std::runtime_error("Tensoriality of the physical quantity not supported.");
        };
    }


private:
    template<typename TensorPhysicalQuantityT>
    struct _Key
    {
        BdryLoc::Type location;
        TensorPhysicalQuantityT quantity;

        // this enable to automatically define comparison operator so that the Key struct can be
        // used as a key in a std::map
        auto operator<=>(_Key const&) const = default;
    };

    using _scalar_key_type = _Key<typename PhysicalQuantityT::Scalar>;
    using _vector_key_type = _Key<typename PhysicalQuantityT::Vector>;

    using _scalar_condition_map_type
        = std::map<_scalar_key_type, std::shared_ptr<scalar_condition_type>>;
    using _vector_condition_map_type
        = std::map<_vector_key_type, std::shared_ptr<vector_condition_type>>;

    _scalar_condition_map_type scalar_conditions_;
    _vector_condition_map_type vector_conditions_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER
