#ifndef PHARE_CORE_NUMERICS_FARADAY_FARADAY_TRANSFORMER_HPP
#define PHARE_CORE_NUMERICS_FARADAY_FARADAY_TRANSFORMER_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/differential_operators.hpp"

#include <array>
#include <stdexcept>

namespace PHARE::core
{

// what one tile's block needs, built on the host: plain views, copyable to the device
template<typename Curl_t, typename VecField_t, typename Box_t>
struct FaradayTile
{
    Curl_t curl;
    VecField_t B, Bnew;
    std::array<Box_t, 3> boxes;
};


// Bnew = B - dt curl(E) over one patch: Faraday per point, on the domain box of each Bnew component
//  untiled: the patch's fields, or one tile's views (on_tile)
//  tiled:   every tile of the patch, CPU tile by tile, GPU one launch with a block per tile
//  Bnew ghosts between tiles are not synced, see sync_inner_ghosts
//  B and Bnew may be the same vecfield (corrector), curl(E) does not read B
class FaradayTransformer
{
public:
    // the box Bnew component Bnc is computed on, as layout.evalOnBox below
    template<typename GridLayout, typename Field>
    static auto box_for(GridLayout const& layout, Field const& Bnc)
    {
        return layout.domainBoxFor(Bnc);
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double const dt) const
        requires(not has_tiled_field_type_c<VecField>)
    {
        auto const curlE = curl(layout, E);

        // captureless kernel: fields travel as arguments so the lambda is device-copyable
        auto const fn = [] _PHARE_ALL_FN_(auto const& ijk, auto& Bnc, auto const& Bc,
                                          auto const& curlEc, double const dt_) {
            Faraday{}(ijk, Bnc, Bc, curlEc, dt_);
        };

        auto const eval = [&](auto& Bnc, auto const& Bc, auto const& curlEc) {
            layout.evalOnBox(Bnc, fn, Bnc, Bc, curlEc, dt);
        };

        eval(Bnew(Component::X), B(Component::X), curlE.x);
        eval(Bnew(Component::Y), B(Component::Y), curlE.y);
        eval(Bnew(Component::Z), B(Component::Z), curlE.z);
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& /*patch layout*/, VecField const& B, VecField const& E,
                    VecField& Bnew, double const dt) const
        requires(has_tiled_field_type_c<VecField>)
    {
        if constexpr (VecField::field_type::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_tiles_gpu(B, E, Bnew, dt);
        else
            for (std::size_t ti = 0; ti < tile_count(Bnew); ++ti)
                on_tile(ti, B, E, Bnew, dt);
    }

    // one tile of tiled B, E and Bnew, eg as one task of a thread pool
    void on_tile(std::size_t const tile_idx, auto const& B, auto const& E, auto& Bnew,
                 double const dt) const
    {
        tile_exec_with_layout_at(
            tile_idx,
            [&](auto& layout, auto const& B_t, auto const& E_t, auto Bnew_t) {
                (*this)(layout, B_t, E_t, Bnew_t, dt);
            },
            B, E, Bnew);
    }

private:
    // block per tile, block threads stride each Bnew component's box
    //  waits for the launch before returning: the per tile items live in this scope
    //   todo: async (no wait), the items must then outlive the launch, see .gpu_todo
    template<typename VecField>
    void on_tiles_gpu([[maybe_unused]] VecField const& B, [[maybe_unused]] VecField const& E,
                      [[maybe_unused]] VecField& Bnew, [[maybe_unused]] double const dt) const
    {
#if PHARE_HAVE_GPU
        using Launcher = gpu::FieldTileLauncher</*sync=*/true>;

        auto const make_item = [&](std::size_t const ti) {
            auto const& layout = Bnew[0][ti].layout();
            auto const Bt      = tile_at(B, ti);
            auto const Et      = tile_at(E, ti);
            auto Bnt           = tile_at(Bnew, ti);
            auto const boxes = for_N_make_array<3>([&](auto c) { return box_for(layout, Bnt[c]); });
            auto const curlE = curl(layout, Et);
            using Box_t      = std::decay_t<decltype(boxes)>::value_type;
            using Item_t     = FaradayTile<std::decay_t<decltype(curlE)>, decltype(Bnt), Box_t>;
            return Item_t{curlE, Bt, Bnt, boxes};
        };
        using Item_t = decltype(make_item(0));

        auto const n_tiles = tile_count(Bnew);
        gpu::Vec_t<Item_t> items;
        items.reserve(n_tiles);
        for (std::size_t ti = 0; ti < n_tiles; ++ti)
            items.emplace_back(make_item(ti));

        auto* const ip = items.data();
        Launcher{n_tiles}([=] _PHARE_ALL_FN_() mutable {
            auto& item = ip[Launcher::tile_idx()];
            Launcher::for_each_point(item.boxes[0], [&](auto const& ijk) {
                Faraday{}(ijk, item.Bnew[0], item.B[0], item.curl.x, dt);
            });
            Launcher::for_each_point(item.boxes[1], [&](auto const& ijk) {
                Faraday{}(ijk, item.Bnew[1], item.B[1], item.curl.y, dt);
            });
            Launcher::for_each_point(item.boxes[2], [&](auto const& ijk) {
                Faraday{}(ijk, item.Bnew[2], item.B[2], item.curl.z, dt);
            });
        });
#else
        throw std::runtime_error("FaradayTransformer: GPU_UNIFIED tiles need a GPU build");
#endif
    }
};

} // namespace PHARE::core

#endif /* PHARE_CORE_NUMERICS_FARADAY_FARADAY_TRANSFORMER_HPP */
