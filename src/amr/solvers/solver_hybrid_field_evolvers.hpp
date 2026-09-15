#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"

#include "amr/resources_manager/amr_utils.hpp"

namespace PHARE::solver
{

class OhmSingleTransformer
{
    using info_type = core::OhmInfo;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

public:
    explicit OhmSingleTransformer(info_type const& info)
        : info_{info}
    {
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew)
    {
        if constexpr (core::is_field_tile_set_v<Field>)
        {
            using Tile_vt = Field::value_type::value_type;
            static_assert(core::is_field_v<Tile_vt>);
            using V_t = core::basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < n().size(); ++tidx)
            {
                auto Enw             = Enew.template as<V_t>([&](auto& c) { return c()[tidx](); });
                auto const& tile_lay = n()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                core::Ohm<TL>{info_, tile_lay}(n()[tidx](), tt<V_t>(Ve, tidx), Pe()[tidx](),
                                               tt<V_t>(B, tidx), tt<V_t>(J, tidx), Enw);
            }
            for (std::uint8_t i = 0; i < 3; ++i)
                Enew[i].sync_inner_ghosts();
        }
        else
        {
            core::Ohm<GridLayout>{info_, layout}(n, Ve, Pe, B, J, Enew);
        }
    }

    info_type info_;
};

template<typename Model>
class OhmLevelTransformer : public OhmSingleTransformer
{
    using Super      = OhmSingleTransformer;
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using info_type  = core::OhmInfo;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

public:
    explicit OhmLevelTransformer(info_type const& info, level_t& level, Model& model)
        : Super{info}
        , level_{level}
        , model_{model}
    {
    }

    void operator()(auto& B, auto& J, auto& E, auto& electrons)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto& n     = electrons.density();
            auto& Ve    = electrons.velocity();
            auto& Pe    = electrons.pressure();
            Super::operator()(layout, n, Ve, Pe, B, J, E);
        }
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    level_t& level_;
    Model& model_;
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;

} // namespace PHARE::solver


#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
