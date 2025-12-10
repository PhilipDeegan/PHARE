#ifndef PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER
#define PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"

#include "initializer/data_provider.hpp"

#include <map>
#include <memory>
#include <initializer_list>

namespace PHARE::core
{

template<IsField FieldT, IsGridLayout GridLayoutT>
class BoundaryManager
{
public:
    using ConditionPtr      = std::shared_ptr<IFieldBoundaryCondition<FieldT, GridLayoutT>>;
    using PhysicalQuantityT = FieldT::physical_quantity_type;

    struct Key
    {
        BdryLoc::Type location;
        PhysicalQuantityT quantity;

        // this enable to automatically define comparison operator so that the Key struct can be
        // used as a key in a std::map
        auto operator<=>(Key const&) const = default;
    };

    using ConditionMap = std::map<Key, ConditionPtr>;

    BoundaryManager() = delete;
    BoundaryManager(PHARE::initializer::PHAREDict const& dict,
                    std::initializer_list<PhysicalQuantityT> const quantities)
    {
        // do not actually read anything, and set everyone to Neumann for the moment
        auto locations = {BdryLoc::Type::XLO, BdryLoc::Type::XHI};
        for (auto const location : locations)
        {
            for (auto const quantity : quantities)
            {
                // we will need a factory here
                boundary_conditions_[Key({location, quantity})]
                    = std::make_shared<FieldNeumannBoundaryCondition<FieldT, GridLayoutT>>(
                        location, quantity);
            };
        };
    }

    ConditionPtr getBoundaryCondition(BdryLoc::Type location, PhysicalQuantityT quantity) const
    {
        auto it = boundary_conditions_.find({location, quantity});
        return (it != boundary_conditions_.end()) ? it->second : nullptr;
    }


private:
    ConditionMap boundary_conditions_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_MHD_BOUNDARY_MANAGER
