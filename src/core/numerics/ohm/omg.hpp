#ifndef PHARE_CORE_NUMERICS_OHM_OMG_HPP
#define PHARE_CORE_NUMERICS_OHM_OMG_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/grid_tiles.hpp"

#include "core/data/tensorfield/tensorfield.hpp"


namespace PHARE::core
{


class OhmSingleTransformer
{
    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i]; });
    }

public:
    explicit OhmSingleTransformer(initializer::PHAREDict const& dict)
        : eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
        , hyper_mode{cppdict::get_value(dict, "hyper_mode", std::string{"constant"}) == "constant"
                         ? HyperMode::constant
                         : HyperMode::spatial}
    {
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew)
    {
        using core_type = Ohm_ref<GridLayout>;

        if constexpr (is_field_tile_set_v<Field>)
        {
            using Tile_vt = Field::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < n().size(); ++tidx)
            {
                auto Enw = Enew.template as<V_t>([&](auto& c) { return c()[tidx]; });
                core_type{n()[tidx].layout(), eta_, nu_, hyper_mode}(n()[tidx](), tt<V_t>(Ve, tidx),
                                                                     Pe()[tidx](), tt<V_t>(B, tidx),
                                                                     tt<V_t>(J, tidx), Enw);
            }
            // ?
            for (std::uint8_t i = 0; i < 3; ++i)
                Enew[i].sync_inner_ghosts();
        }
        else
        {
            core_type{layout, eta_, nu_, hyper_mode}(n, Ve, Pe, B, J, Enew);
        }
    }

    double const eta_;
    double const nu_;
    HyperMode const hyper_mode;
};


} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_OHM_OMG_HPP