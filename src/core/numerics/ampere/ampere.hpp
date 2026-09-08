#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP


#include "core/data/grid/grid_tiles.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/types.hpp"


namespace PHARE::core
{



template<typename GridLayout>
class Ampere
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Ampere(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField& J) _PHARE_ALL_FN_
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        Point<std::uint32_t, dimension> const shrink{ConstArray<std::uint32_t, dimension>(1)};

        layout_.evalOnShrinkedGhostBox(
            Jx, shrink, [] _PHARE_ALL_FN_(auto&&... args) { JxEq_(args...); }, Jx, B, layout_);
        layout_.evalOnShrinkedGhostBox(
            Jy, shrink, [] _PHARE_ALL_FN_(auto&&... args) { JyEq_(args...); }, Jy, B, layout_);
        layout_.evalOnShrinkedGhostBox(
            Jz, shrink, [] _PHARE_ALL_FN_(auto&&... args) { JzEq_(args...); }, Jz, B, layout_);
    }


private:
    GridLayout layout_;


    template<typename IJK, typename... Args>
    static void JxEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto&& [Jx, B, layout] = std::forward_as_tuple(args...);
        auto&& [_, By, Bz]     = B();

        if constexpr (dimension == 1)
            Jx(ijk) = 0.0;

        if constexpr (dimension == 2)

            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, ijk);

        if constexpr (dimension == 3)
            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, ijk)
                      - layout.template deriv<Direction::Z>(By, ijk);
    }

    template<typename IJK, typename... Args>
    static void JyEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto&& [Jy, B, layout] = std::forward_as_tuple(args...);
        auto&& [Bx, By, Bz]    = B();

        if constexpr (dimension == 1 || dimension == 2)
            Jy(ijk) = -layout.template deriv<Direction::X>(Bz, ijk);

        if constexpr (dimension == 3)
            Jy(ijk) = layout.template deriv<Direction::Z>(Bx, ijk)
                      - layout.template deriv<Direction::X>(Bz, ijk);
    }

    template<typename IJK, typename... Args>
    static void JzEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto&& [Jz, B, layout] = std::forward_as_tuple(args...);
        auto&& [Bx, By, Bz]    = B();

        if constexpr (dimension == 1)
            Jz(ijk) = layout.template deriv<Direction::X>(By, ijk);

        else
            Jz(ijk) = layout.template deriv<Direction::X>(By, ijk)
                      - layout.template deriv<Direction::Y>(Bx, ijk);
    }
};


class AmpereSingleTransformer
{
    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i]; });
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
    {
        using field_type = VecField::field_type;

        if constexpr (is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < J[0]().size(); ++tidx)
            {
                auto Jt              = J.template as<V_t>([&](auto& c) { return c()[tidx]; });
                auto const& tile_lay = J[0]()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                Ampere<TL>{tile_lay}(tt<V_t>(B, tidx), Jt);
            }
            for (std::uint8_t i = 0; i < 3; ++i)
                J[i].sync_inner_ghosts();
        }
        else
        {
            Ampere<GridLayout>{layout}(B, J);
        }

        check_tensor_field(J, layout);
    }
};

} // namespace PHARE::core
#endif
