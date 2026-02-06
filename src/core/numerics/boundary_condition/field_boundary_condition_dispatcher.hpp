#ifndef PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER
#define PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"


namespace PHARE::core
{
/**
 * @brief Intermediate dispatcher for field boundary conditions class, inheriting from @link
 * PHARE::core::IFieldBoundaryCondition @endlink, and whom concrete implementations must inherit
 * from.
 *
 * Provides a mechanism to dispatch runtime boundary information (location, centering)
 * to compile-time specialized methods in concrete implementations. It implements the Curious
 * Recurring Template Pattern so the complicated dispatching code is not duplicated in concrete
 * implementations of field boundary conditions. Actual implementations are expected to implement
 * an @c apply_specialized templated function with the following interface:
 *
 * @code
 * template<Direction dir, Side side, QtyCentering... centerings>
 * void apply_specialized(
 *      ScalarOrTensorFieldT& field,
 *      Box<std::uint32_t, dimension> const& local_ghost_box,
 *      GridLayoutT const& grid_layout,
 *      double const& time
 * );
 * @endcode
 *
 * @tparam ScalarOrTensorFieldT Type of field managed.
 * @tparam GridLayoutT Grid layout configuration.
 * @tparam Derived The concrete class inheriting from this dispatcher.
 *
 */
template<typename ScalarOrTensorFieldT, IsGridLayout GridLayoutT, typename Derived>
class FieldBoundaryConditionDispatcher
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using typename Super::field_type;
    using typename Super::physical_quantity_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr bool is_scalar   = Super::is_scalar;
    static constexpr size_t N         = Super::N;

    /**
     * @brief Construct the dispatcher and initialize centering information.
     * @param location Boundary location.
     * @param physical_quantity Physical quantity.
     */
    FieldBoundaryConditionDispatcher(BdryLoc::Type const& location,
                                     physical_quantity_type const& physicalQuantity)
        : location_{location}
        , direction_{BdryLoc::direction(location)}
        , side_{BdryLoc::side(location)}
        , physicalQuantity_{physicalQuantity}
    {
        if constexpr (is_scalar)
        {
            centerings_[0]
                = GridLayoutT::centering(physicalQuantity)[static_cast<size_t>(direction_)];
        }
        else
        {
            auto full_centerings = GridLayoutT::centering(physicalQuantity);
            for (size_t i = 0; i < N; ++i)
                centerings_[i] = full_centerings[i][static_cast<size_t>(direction_)];
        }
    };

    BdryLoc::Type getLocation() const override { return location_; };

    Direction getDirection() const override { return direction_; };

    Side getSide() const override { return side_; };

    physical_quantity_type getPhysicalQuantity() const override { return physicalQuantity_; };

    std::array<QtyCentering, N> getQtyCenterings() const override { return centerings_; };

    /**
     * @brief Implements the @link PHARE::core::IFieldBoundaryCondition::apply @endlink abstract
     * function.
     *
     * Triggers the recursive dispatching of centerings, directions, and sides to
     * specialized implementations.
     */
    void apply(ScalarOrTensorFieldT& field, Box<std::uint32_t, dimension> const& localGhostBox,
               GridLayoutT const& gridLayout, double const& time) override
    {
        dispatch_centerings<>(field, localGhostBox, gridLayout, time);
    }

protected:
    /**
     * @brief Helper that recursively promote runtime centerings, direction and side to compile-time
     * tags.
     *
     * The recursive character of this helper is necessary because we need to promote the centering
     * of all components of the (tensor) field in the direction normal to the boundary.
     */
    template<QtyCentering... AlreadyPromoted>
    void dispatch_centerings(ScalarOrTensorFieldT& field,
                             Box<std::uint32_t, dimension> const& localGhostBox,
                             GridLayoutT const& gridLayout, double const& time)
    {
        constexpr size_t nAlreadyPromoted = sizeof...(AlreadyPromoted);

        if constexpr (nAlreadyPromoted == N)
        // base case: all directional centerings have been promoted
        {
            auto d_v = promote<Direction::X, Direction::Y, Direction::Z>(direction_);
            auto s_v = promote<Side::LOWER, Side::UPPER>(side_);

            std::visit(
                [&](auto d_tag, auto s_tag) {
                    static_cast<Derived*>(this)
                        ->template apply_specialized<d_tag.value, s_tag.value, AlreadyPromoted...>(
                            field, localGhostBox, gridLayout, time);
                },
                d_v, s_v);
        }
        else
        // we grow the list of promoted centerings to call the next recursive version of the
        // function
        {
            auto c_v
                = promote<QtyCentering::primal, QtyCentering::dual>(centerings_[nAlreadyPromoted]);

            std::visit(
                [&](auto c_tag) {
                    this->dispatch_centerings<AlreadyPromoted..., c_tag.value>(field, localGhostBox,
                                                                               gridLayout, time);
                },
                c_v);
        };
    }


private:
    BdryLoc::Type location_;
    Direction direction_;
    Side side_;
    physical_quantity_type physicalQuantity_;
    std::array<QtyCentering, N> centerings_;
};
} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DISPATCHER
