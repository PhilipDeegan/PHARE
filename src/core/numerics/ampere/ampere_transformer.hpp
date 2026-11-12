#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_TRANSFORMER_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_TRANSFORMER_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/kernels/region_chunks.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/differential_operators.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

namespace PHARE::core
{

// what one tile's block needs, built on the host: plain views, copyable to the device
template<typename Curl_t, typename VecField_t, typename Box_t>
struct AmpereTile
{
    Curl_t curl;
    VecField_t J;
    std::array<Box_t, 3> boxes;
};


// J = curl(B) over one patch: Ampere per point, on the box `box` selects per J component
//  untiled: the patch's fields, or one tile's views (on_tile)
//  tiled:   every tile of the patch, CPU tile by tile, GPU one launch with a block per tile
//  J ghosts between tiles are not synced, see sync_inner_ghosts
template<AmpereBox box = AmpereBox{}>
class AmpereTransformer
{
public:
    // the box J component Jc is computed on
    //  must match the layout evaluator used in the untiled operator() below
    template<typename GridLayout, typename Field>
    static auto box_for(GridLayout const& layout, Field const& Jc)
    {
        if constexpr (box.mode == AmpereMode::ShrinkedGhost)
            return shrink(layout.ghostBoxFor(Jc), amount<GridLayout>());
        else
            return grow(layout.domainBoxFor(Jc), amount<GridLayout>());
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J) const
        requires(not has_tiled_field_type_c<VecField>)
    {
        if constexpr (VecField::field_type::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            auto items = gpu_items(layout, B, J);
            run_gpu(items);
            return;
        }

        auto const curlB = curl(layout, B);

        // captureless kernel: fields travel as arguments so the lambda is device-copyable
        auto const fn = [] _PHARE_ALL_FN_(auto const& ijk, auto& Jc, auto const& curlc) {
            Ampere{}(ijk, Jc, curlc);
        };

        auto const eval = [&](auto& Jc, auto const& curlc) {
            if constexpr (box.mode == AmpereMode::ShrinkedGhost)
                layout.evalOnShrinkedGhostBox(Jc, amount<GridLayout>(), fn, Jc, curlc);
            else
                layout.evalOnBiggerBox(Jc, amount<GridLayout>(), fn, Jc, curlc);
        };

        eval(J(Component::X), curlB.x);
        eval(J(Component::Y), curlB.y);
        eval(J(Component::Z), curlB.z);
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& /*patch layout*/, VecField const& B, VecField& J) const
        requires(has_tiled_field_type_c<VecField>)
    {
        if constexpr (VecField::field_type::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_tiles_gpu(B, J);
        else
            for (std::size_t ti = 0; ti < tile_count(J); ++ti)
                on_tile(ti, B, J);
    }

    // one tile of tiled B and J, eg as one task of a thread pool
    void on_tile(std::size_t const tile_idx, auto const& B, auto& J) const
    {
        tile_exec_with_layout_at(
            tile_idx, [&](auto& layout, auto const& B_t, auto J_t) { (*this)(layout, B_t, J_t); },
            B, J);
    }

private:
    template<typename GridLayout>
    static auto amount()
    {
        auto constexpr dim = GridLayout::dimension;
        return Point<std::uint32_t, dim>{ConstArray<std::uint32_t, dim>(box.amount)};
    }

public:
    // GPU_UNIFIED, as region items (an untiled patch is one region, a tiled one a region per tile)
    //  collected per patch or across a level (see AmpereLevelTransformer), one launch per run_gpu
    template<typename GridLayout, typename VecField>
    static auto region_layout(GridLayout const& layout, VecField const& J, std::size_t const ti)
    {
        if constexpr (has_tiled_field_type_c<VecField>)
            return J[0][ti].layout();
        else
            return layout;
    }

    template<typename VecField>
    static auto region_view(VecField&& vf, std::size_t const ti)
    {
        if constexpr (has_tiled_field_type_c<std::remove_cvref_t<VecField>>)
            return tile_at(vf, ti);
        else
            return untiled_view(vf);
    }

    template<typename GridLayout, typename B_t, typename J_t>
    static auto make_item(GridLayout const& layout, B_t const& B, J_t const& J)
    {
        auto const boxes = for_N_make_array<3>([&](auto c) { return box_for(layout, J[c]); });
        auto const curlB = curl(layout, B);
        using Box_t      = std::decay_t<decltype(boxes)>::value_type;
        return AmpereTile<std::decay_t<decltype(curlB)>, J_t, Box_t>{curlB, J, boxes};
    }

    template<typename GridLayout, typename VecField>
    using item_t = decltype(make_item(
        std::declval<GridLayout const&>(), region_view(std::declval<VecField const&>(), 0),
        region_view(std::declval<VecField&>(), 0)));

    // appends one patch's items
    template<typename GridLayout, typename VecField>
    static void gpu_items(GridLayout const& layout, VecField const& B, VecField& J, auto& items)
    {
        auto const n_regions = [&]() -> std::size_t {
            if constexpr (has_tiled_field_type_c<VecField>)
                return tile_count(J);
            else
                return 1;
        }();
        for (std::size_t ti = 0; ti < n_regions; ++ti)
            items.emplace_back(make_item(region_layout(layout, J, ti), region_view(B, ti),
                                         region_view(J, ti)));
    }

    template<typename GridLayout, typename VecField>
    static auto gpu_items(GridLayout const& layout, VecField const& B, VecField& J)
    {
        gpu_item_list_t<GridLayout, VecField> items;
        gpu_items(layout, B, J, items);
        return items;
    }

#if PHARE_HAVE_GPU
    template<typename GridLayout, typename VecField>
    using gpu_item_list_t = gpu::Vec_t<item_t<GridLayout, VecField>>;
#else
    template<typename GridLayout, typename VecField>
    using gpu_item_list_t = std::vector<item_t<GridLayout, VecField>>;
#endif

    static void run_gpu(auto& items)
    {
        run_region_chunks(items, [] _PHARE_ALL_FN_(auto& item, auto const c, auto const& ijk) {
            if (c == 0)
                Ampere{}(ijk, item.J[0], item.curl.x);
            else if (c == 1)
                Ampere{}(ijk, item.J[1], item.curl.y);
            else
                Ampere{}(ijk, item.J[2], item.curl.z);
        });
    }

    template<typename VecField>
    void on_tiles_gpu(VecField const& B, VecField& J) const
    {
        auto items = gpu_items(J[0][0].layout(), B, J);
        run_gpu(items);
    }
};

} // namespace PHARE::core

#endif /* PHARE_CORE_NUMERICS_AMPERE_AMPERE_TRANSFORMER_HPP */
