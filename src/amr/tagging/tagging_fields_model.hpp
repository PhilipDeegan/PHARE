#ifndef PHARE_TAGGING_FIELDS_MODEL_HPP
#define PHARE_TAGGING_FIELDS_MODEL_HPP

#include "core/def.hpp"
#include "phare_core.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "initializer/data_provider.hpp"

#include "amr/tagging/tagger_strategy.hpp"
#include "amr/tagging/default_tagger_strategy.hpp"
#include "amr/tagging/gradient_scale_tagger_strategy.hpp"
#include "amr/tagging/lohner_per_direction_tagger_strategy.hpp"
#include "amr/tagging/lohner_cross_derivatives_tagger_strategy.hpp"

#include <array>
#include <memory>
#include <string>
#include <span>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace PHARE::amr::tagging
{
/**
 * @brief A hybrid model made of its tagged fields only, viewing buffers it does not own: what the
 * tagging strategies read of a model (its B, its layout type), so they can tag fields computed
 * outside a simulation, e.g. a run's initial fields evaluated from its input functions.
 *
 * model_type_name makes solver::is_hybrid_model_v true for it: the strategies resolve its
 * quantities as a hybrid model's (state.electromag.B).
 *
 * @tparam opts the simulation options (dimension, interpolation order) giving the layout
 */
template<auto opts>
struct HybridTaggingFields
{
    static constexpr std::string_view model_type_name = "HybridModel";
    static constexpr auto dimension                   = opts.dimension;

    using types           = core::PHARE_Types<opts>::Hybrid;
    using gridlayout_type = types::GridLayout_t;
    using field_type      = types::Field_t;
    using vecfield_type   = types::VecField_t;

    struct
    {
        struct
        {
            vecfield_type B{"EM_B", core::HybridQuantity::Vector::B};
        } electromag;
    } state;
};


/**
 * @brief The tagging strategy of a run's tagging config, for a model type.
 *
 * @param dict the config at "simulation/AMR/refinement/tagging", with "hybrid_method" for a
 * hybrid model (as TaggerFactory reads it)
 */
template<typename Model>
std::unique_ptr<TaggerStrategy<Model>> makeTaggerStrategy(initializer::PHAREDict const& dict)
{
    std::string const method = dict["hybrid_method"].template to<std::string>();
    if (method == "default")
        return std::make_unique<DefaultTaggerStrategy<Model>>(dict);
    if (method == "lohner-per-direction")
        return std::make_unique<LohnerPerDirectionTaggerStrategy<Model>>(dict);
    if (method == "lohner-cross-derivatives")
        return std::make_unique<LohnerCrossDerivativesTaggerStrategy<Model>>(dict);
    if (method == "gradient-scale")
        return std::make_unique<GradientScaleTaggerStrategy<Model>>(dict);
    throw std::runtime_error("tagging: unknown hybrid_method '" + method + "'");
}


/**
 * @brief The tags of one patch from its B, outside a simulation.
 *
 * @param dict the tagging config, as for makeTaggerStrategy
 * @param meshSize the patch level's mesh size
 * @param lower, upper the patch's cells (AMR index space of its level)
 * @param B the components' values over their allocation (ghosts included, the layout's
 * allocSize, C ordering)
 * @return one tag per cell of the patch, Fortran ordering (as SAMRAI's tag arrays)
 */
template<auto opts>
std::vector<int>
tagPatch(initializer::PHAREDict const& dict, std::array<double, opts.dimension> const& meshSize,
         std::array<int, opts.dimension> const& lower, std::array<int, opts.dimension> const& upper,
         std::array<std::span<double>, 3> const& B)
{
    using Model          = HybridTaggingFields<opts>;
    auto constexpr dim   = opts.dimension;
    using Scalar         = core::HybridQuantity::Scalar;
    using field_type     = Model::field_type;
    using layout_type    = Model::gridlayout_type;
    auto constexpr names = std::array{"EM_B_x", "EM_B_y", "EM_B_z"};
    auto constexpr qties = std::array{Scalar::Bx, Scalar::By, Scalar::Bz};

    std::array<std::uint32_t, dim> cells;
    core::Point<double, dim> origin;
    for (std::size_t i = 0; i < dim; ++i)
    {
        if (upper[i] < lower[i])
            throw std::runtime_error("tagging: empty patch");
        cells[i]  = static_cast<std::uint32_t>(upper[i] - lower[i] + 1);
        origin[i] = lower[i] * meshSize[i];
    }
    layout_type const layout{
        meshSize, cells, origin,
        core::Box<int, dim>{core::Point<int, dim>{lower}, core::Point<int, dim>{upper}}};

    std::vector<field_type> fields;
    fields.reserve(3);
    for (std::size_t c = 0; c < 3; ++c)
    {
        auto const shape = layout.allocSize(qties[c]);
        std::size_t size = 1;
        for (auto const s : shape)
            size *= s;
        if (B[c].size() != size)
            throw std::runtime_error(std::string{"tagging: "} + names[c] + " has "
                                     + std::to_string(B[c].size()) + " values, its layout "
                                     + std::to_string(size));
        fields.emplace_back(names[c], qties[c], B[c].data(), shape);
    }

    Model model;
    model.state.electromag.B.setBuffer(&fields);

    std::size_t n = 1;
    for (auto const c : cells)
        n *= c;
    std::vector<int> tags(n, 0);
    makeTaggerStrategy<Model>(dict)->tag(model, layout, tags.data());
    return tags;
}

} // namespace PHARE::amr::tagging

#endif
