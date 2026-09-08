#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP


#include "core/def.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


namespace PHARE::core
{


template<typename GridLayout>
class Faraday
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Faraday(GridLayout const& layout)
        : layout_{layout}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt) _PHARE_ALL_FN_
    {
        dt_ = dt;
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        layout_.evalOnBox(
            Bxnew, [] _PHARE_ALL_FN_(auto&&... args) { BxEq_(args...); }, Bx, E, Bxnew, *this);
        layout_.evalOnBox(
            Bynew, [] _PHARE_ALL_FN_(auto&&... args) { ByEq_(args...); }, By, E, Bynew, *this);
        layout_.evalOnBox(
            Bznew, [] _PHARE_ALL_FN_(auto&&... args) { BzEq_(args...); }, Bz, E, Bznew, *this);
    }

private:
    GridLayout layout_;
    double dt_;



    template<typename IJK, typename... Args>
    static void BxEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bx, E, Bxnew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [_, Ey, Ez]          = E();

        if constexpr (dimension == 1)
            Bxnew(ijk) = Bx(ijk);

        if constexpr (dimension == 2)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk);

        if constexpr (dimension == 3)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk)
                         + dt * layout.template deriv<Direction::Z>(Ey, ijk);
    }

    template<typename IJK, typename... Args>
    static void ByEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [By, E, Bynew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, _, Ez]          = E();

        if constexpr (dimension == 1 || dimension == 2)
            Bynew(ijk) = By(ijk) + dt * layout.template deriv<Direction::X>(Ez, ijk);

        if constexpr (dimension == 3)
            Bynew(ijk) = By(ijk) - dt * layout.template deriv<Direction::Z>(Ex, ijk)
                         + dt * layout.template deriv<Direction::X>(Ez, ijk);
    }

    template<typename IJK, typename... Args>
    static void BzEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bz, E, Bznew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, Ey, _]          = E();

        if constexpr (dimension == 1)
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk);

        else
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk)
                         + dt * layout.template deriv<Direction::Y>(Ex, ijk);
    }
};


class FaradaySingleTransformer
{
    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E,
                    VecField& Bnew, double dt)
    {
        using field_type = VecField::field_type;

        if constexpr (is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < B[0]().size(); ++tidx)
            {
                auto Bnw             = Bnew.template as<V_t>([&](auto& c) { return c()[tidx](); });
                auto const& tile_lay = B[0]()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                Faraday<TL>{tile_lay}(tt<V_t>(B, tidx), tt<V_t>(E, tidx), Bnw, dt);
            }
            for (std::uint8_t i = 0; i < 3; ++i)
                Bnew[i].sync_inner_ghosts();
        }
        else
        {
            Faraday<GridLayout>{layout}(B, E, Bnew, dt);
        }
    }
};

} // namespace PHARE::core


#endif
