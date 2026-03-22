#ifndef PHARE_CORE_NUMERICS_AMPERE_AMDERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMDERE_HPP


#include "core/data/grid/grid_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"

#include "core/data/tensorfield/tensorfield.hpp"



namespace PHARE::core
{



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
        using core_type  = Ampere_ref<GridLayout>;
        using field_type = VecField::field_type;


        if constexpr (is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < J[0]().size(); ++tidx)
            {
                auto Jt = J.template as<V_t>([&](auto& c) { return c()[tidx]; });
                core_type{J[0]()[tidx].layout()}(tt<V_t>(B, tidx), Jt);
            }

            for (std::uint8_t i = 0; i < 3; ++i)
                J[i].sync_inner_ghosts();
        }
        else
        {
            core_type{layout}(B, J);
        }
    }
};



} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_AMPERE_AMDERE_HPP