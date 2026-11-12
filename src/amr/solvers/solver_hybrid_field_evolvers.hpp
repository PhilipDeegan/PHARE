#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm_transformer.hpp"
#include "core/data/field/field.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/data/field/field_tiles_gpu.hpp"

#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/resources_manager/amr_utils.hpp"

namespace PHARE::solver
{

template<typename Model>
class OhmLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using info_type  = core::OhmInfo;
    using core_type  = core::OhmTransformer;

public:
    explicit OhmLevelTransformer(info_type const& info, level_t& level, Model& model)
        : ohm_{info}
        , level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(VecField& B, VecField& J, VecField& E, auto& electrons)
        requires(not core::has_tiled_field_type_c<VecField>)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto& n     = electrons.density();
            auto& Ve    = electrons.velocity();
            auto& Pe    = electrons.pressure();
            ohm_(layout, n, Ve, Pe, B, J, E);
        }
    }

    template<typename VecField>
    void operator()(VecField& B, VecField& J, VecField& E, auto& electrons)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        auto accessor     = amr::make_model_level_accessor(level_, model_, electrons, B, J, E);
        auto const output = [](auto& view) -> auto& { return std::get<3>(view.args); };

        if constexpr (VecField::field_type::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            // per patch: one launch, a block per tile, then the tile ghost sync on device
            amr::patch_level_exec(accessor, [&](auto& view) {
                auto& [electrons_v, B_v, J_v, E_v] = view.args;
                ohm_(view.layout, electrons_v.density(), electrons_v.velocity(),
                     electrons_v.pressure(), B_v, J_v, E_v);
                core::sync_inner_ghosts_on_device(output(view));
            });
        }
        else
            tiled_level_transform(accessor, output, [this](auto& view, auto const tile_idx) {
                auto& [electrons_v, B_v, J_v, E_v] = view.args;
                ohm_.on_tile(tile_idx, electrons_v.density(), electrons_v.velocity(),
                             electrons_v.pressure(), B_v, J_v, E_v);
            });
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    core_type const ohm_;
    level_t& level_;
    Model& model_;
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;


// Ve and Pe from the ions moments and J, per patch or per tile on the thread pools
template<typename Model>
class ElectronsLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using field_type = Model::field_type;

public:
    explicit ElectronsLevelTransformer(level_t& level, Model& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(auto& electrons)
        requires(not core::is_field_tile_set_v<field_type>)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons))
            electrons.update(amr::layoutFromPatch<GridLayout>(*patch));
    }

    void operator()(auto& electrons)
        requires(core::is_field_tile_set_v<field_type>)
    {
        auto accessor = amr::make_model_level_accessor(level_, model_, electrons);
        tiled_level_exec(
            accessor, [](auto& view) -> auto& { return std::get<0>(view.args).density(); },
            [](auto& view, auto const tile_idx) { std::get<0>(view.args).update(tile_idx); });
    }

    level_t& level_;
    Model& model_;
};

template<typename Model>
ElectronsLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> ElectronsLevelTransformer<Model>;

} // namespace PHARE::solver

#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
