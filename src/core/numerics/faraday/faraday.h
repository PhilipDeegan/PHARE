#ifndef PHARE_FARADAY_H
#define PHARE_FARADAY_H

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"


namespace PHARE::core
{
template<typename GridLayout>
struct StandardFaradayComputer
{
    constexpr static auto dimension = GridLayout::dimension;

    template<typename VecField, typename Field, typename... Indexes>
    void bx(Field const& Bx, VecField const& E, Field& Bxnew,
            Indexes const&... ijk) const _PHARE_ALL_FN_
    {
        auto const& [_, Ey, Ez] = E();

        if constexpr (dimension == 1)
            Bxnew(ijk...) = Bx(ijk...);

        if constexpr (dimension == 2)
            Bxnew(ijk...)
                = Bx(ijk...)
                  - dt * layout.deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::Y>{});

        if constexpr (dimension == 3)
            Bxnew(ijk...) = Bx(ijk...)
                            - dt * layout.deriv(Ez, {ijk...}, DirectionTag<Direction::Y>{})
                            + dt * layout.deriv(Ey, {ijk...}, DirectionTag<Direction::Z>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void by(Field const& By, VecField const& E, Field& Bynew,
            Indexes const&... ijk) const _PHARE_ALL_FN_
    {
        auto const& [Ex, _, Ez] = E();

        if constexpr (dimension == 1)
            Bynew(ijk...)
                = By(ijk...)
                  + dt * layout.deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 2)
            Bynew(ijk...)
                = By(ijk...)
                  + dt * layout.deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 3)
            Bynew(ijk...) = By(ijk...)
                            - dt * layout.deriv(Ex, {ijk...}, DirectionTag<Direction::Z>{})
                            + dt * layout.deriv(Ez, {ijk...}, DirectionTag<Direction::X>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void bz(Field const& Bz, VecField const& E, Field& Bznew,
            Indexes const&... ijk) const _PHARE_ALL_FN_
    {
        auto const& [Ex, Ey, _] = E();

        if constexpr (dimension == 1)
            Bznew(ijk...)
                = Bz(ijk...) - dt * layout.deriv(Ey, {ijk...}, DirectionTag<Direction::X>{});

        else
            Bznew(ijk...) = Bz(ijk...)
                            - dt * layout.deriv(Ey, {ijk...}, DirectionTag<Direction::X>{})
                            + dt * layout.deriv(Ex, {ijk...}, DirectionTag<Direction::Y>{});
    }

    double dt;
    GridLayout& layout;
};



template<typename GridLayout, typename Computer = StandardFaradayComputer<GridLayout>>
class Faraday : public LayoutHolder<GridLayout>
{
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt) _PHARE_ALL_FN_
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");

        if (!(B.isUsable() && E.isUsable() && Bnew.isUsable()))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        Computer op{dt, *this->layout_};
        layout_->scan(Bnew(Component::X), [&]_PHARE_ALL_FN_(auto const&... args) {
            op.bx(B(Component::X), E, Bnew(Component::X), args...);
        });
        layout_->scan(Bnew(Component::Y), [&]_PHARE_ALL_FN_(auto const&... args) {
            op.by(B(Component::Y), E, Bnew(Component::Y), args...);
        });
        layout_->scan(Bnew(Component::Z), [&]_PHARE_ALL_FN_(auto const&... args) {
            op.bz(B(Component::Z), E, Bnew(Component::Z), args...);
        });
    }
};

} // namespace PHARE::core


#endif
