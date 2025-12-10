#ifndef PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER
#define PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"


namespace PHARE::core
{

template<IsField FieldT, IsGridLayout GridLayoutT, typename Derived>
class FieldBoundaryConditionDispatcher : public IFieldBoundaryCondition<FieldT, GridLayoutT>
{
public:
    using Super = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using Super::dimension;
    using typename Super::PhysicalQuantityT;

    FieldBoundaryConditionDispatcher(BdryLoc::Type const& location,
                                   PhysicalQuantityT const& physicalQuantity)
        : location_{location}
        , direction_{BdryLoc::direction(location)}
        , side_{BdryLoc::side(location)}
        , physicalQuantity_{physicalQuantity}
        , centering_{GridLayoutT::centering(
              physicalQuantity)[static_cast<std::size_t>(BdryLoc::direction(location))]}
    {
    } //

    void apply(FieldT& field, Box<std::uint32_t, dimension> const& localGhostBox,
               GridLayoutT const& gridLayout, double const& time) override
    {
        // 1. Lift runtime enums to Variants of Tags
        auto d_v = promote<Direction::X, Direction::Y, Direction::Z>(direction_);
        auto s_v = promote<Side::LOWER, Side::UPPER>(side_);
        auto c_v = promote<QtyCentering::primal, QtyCentering::dual>(centering_);

        // 2. Static Dispatch
        std::visit(
            [&](auto d_tag, auto s_tag, auto c_tag) {
                // Cast 'this' to the Child type to access the template kernel
                static_cast<Derived*>(this)
                    ->template apply_specialized<d_tag.value, s_tag.value, c_tag.value>(
                        field, localGhostBox, gridLayout, time);
            },
            d_v, s_v, c_v);
    }

    BdryLoc::Type getLocation() const override { return location_; };
    Direction getDirection() const override { return direction_; };
    Side getSide() const override { return side_; };
    PhysicalQuantityT getPhysicalQuantity() const override { return physicalQuantity_; };
    QtyCentering getQtyCentering() const override { return centering_; };

private:
    BdryLoc::Type location_;
    Direction direction_;
    Side side_;
    PhysicalQuantityT physicalQuantity_;
    QtyCentering centering_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER
