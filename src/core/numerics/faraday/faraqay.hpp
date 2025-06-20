#ifndef PHARE_FARAQAY_HPP
#define PHARE_FARAQAY_HPP


#include "core/def.hpp"
// #include "core/data/grid/gridlayoutdefs.hpp"
// #include "core/data/grid/gridlayout_utils.hpp"
// #include "core/data/vecfield/vecfield_component.hpp"


#include "core/data/grid/grid_tiles.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "core/data/tensorfield/tensorfield.hpp"



namespace PHARE::core
{



class FaradaySingleTransformer
{
    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i]; });
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double const dt_)
    {
        using core_type  = Faraday_ref<GridLayout>;
        using field_type = VecField::field_type;

        if constexpr (is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < B[0]().size(); ++tidx)
            {
                auto Bnw = Bnew.template as<V_t>([&](auto& c) { return c()[tidx]; });
                core_type{B[0]()[tidx].layout(), dt_}(tt<V_t>(B, tidx), tt<V_t>(E, tidx), Bnw);
            }
        }
        else
        {
            core_type{layout, dt_}(B, E, Bnew);
        }
    }
};



} // namespace PHARE::core

#endif //  PHARE_FARAQAY_HPP
