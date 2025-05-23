#ifndef PHARE_CORE_NUMERICS_OHM_OMG_HPP
#define PHARE_CORE_NUMERICS_OHM_OMG_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/grid_tiles.hpp"


#include <type_traits>

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
    {
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew)
    {
        using core_type = Ohm_ref<GridLayout>;
        // using field_type = std::remove_pointer_t<typename VecField::value_type>::field_type;

        if constexpr (is_field_tile_set_v<Field>)
        {
            // using Vecfield_t    = Electromag_t::vecfield_type;
            // using Field_t       = Vecfield_t::field_type;
            using Tile_vt = Field::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < n().size(); ++tidx)
            {
                auto Enw = Enew.template as<V_t>([&](auto& c) { return c()[tidx]; });
                core_type{n()[tidx].layout(), eta_, nu_}(n()[tidx](), tt<V_t>(Ve, tidx),
                                                         Pe()[tidx](), tt<V_t>(B, tidx),
                                                         tt<V_t>(J, tidx), Enw);
            }
        }
        else
        {
            core_type{layout, eta_, nu_}(n, Ve, Pe, B, J, Enew);
        }
    }

    double const eta_;
    double const nu_;
};


} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_OHM_OMG_HPP