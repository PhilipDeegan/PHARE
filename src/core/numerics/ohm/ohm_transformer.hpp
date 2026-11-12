#ifndef PHARE_CORE_NUMERICS_OHM_OHM_TRANSFORMER_HPP
#define PHARE_CORE_NUMERICS_OHM_OHM_TRANSFORMER_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/differential_operators.hpp"

#include <array>
#include <stdexcept>

namespace PHARE::core
{

// what one tile's block needs, built on the host: plain views, copyable to the device
template<typename Ohm_t, typename Grad_t, typename Field_t, typename VecField_t, typename Box_t>
struct OhmTile
{
    Ohm_t ohm; // holds the tile's layout
    Grad_t gradPe;
    Field_t n;
    VecField_t Ve, B, J, E;
    std::array<Box_t, 3> boxes;
};


// E from Ohm's law over one patch: Ohm per point, on the domain box of each E component
//  untiled: the patch's fields, or one tile's views (on_tile)
//  tiled:   every tile of the patch, CPU tile by tile, GPU one launch with a block per tile
//  E ghosts between tiles are not synced, see sync_inner_ghosts
class OhmTransformer
{
public:
    explicit OhmTransformer(OhmInfo const& info)
        : info_{info}
    {
    }

    // the box E component Ec is computed on, as layout.evalOnBox below
    template<typename GridLayout, typename Field>
    static auto box_for(GridLayout const& layout, Field const& Ec)
    {
        return layout.domainBoxFor(Ec);
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew) const
        requires(not is_field_tile_set_v<Field>)
    {
        with_flags([&]<bool isResistive, bool isHyperResistive>() {
            auto const gradPe = grad(layout, Pe);
            Ohm<GridLayout> const ohm{info_, layout};

            auto const eval = [&]<auto c>(auto const& gradPe_c) {
                // captureless kernel: fields and ohm travel as arguments, device-copyable
                layout.evalOnBox(
                    Enew(c),
                    [] _PHARE_ALL_FN_(auto const& ijk, auto const& ohm_, auto&&... args) {
                        ohm_.template operator()<c, isResistive, isHyperResistive>(ijk, args...);
                    },
                    ohm, n, Ve, gradPe_c, B, J, Enew);
            };

            eval.template operator()<Component::X>(gradPe.x);
            eval.template operator()<Component::Y>(gradPe.y);
            eval.template operator()<Component::Z>(gradPe.z);
        });
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& /*patch layout*/, Field const& n, VecField const& Ve,
                    Field const& Pe, VecField const& B, VecField const& J, VecField& Enew) const
        requires(is_field_tile_set_v<Field>)
    {
        if constexpr (Field::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_tiles_gpu(n, Ve, Pe, B, J, Enew);
        else
            for (std::size_t ti = 0; ti < tile_count(Enew); ++ti)
                on_tile(ti, n, Ve, Pe, B, J, Enew);
    }

    // one tile of the tiled fields, eg as one task of a thread pool
    void on_tile(std::size_t const tile_idx, auto const& n, auto const& Ve, auto const& Pe,
                 auto const& B, auto const& J, auto& Enew) const
    {
        tile_exec_with_layout_at(
            tile_idx,
            [&](auto& layout, auto const& n_t, auto const& Ve_t, auto const& Pe_t, auto const& B_t,
                auto const& J_t, auto E_t) { (*this)(layout, n_t, Ve_t, Pe_t, B_t, J_t, E_t); },
            n, Ve, Pe, B, J, Enew);
    }

private:
    // fn<isResistive, isHyperResistive>(): runtime flags as compile-time tags
    void with_flags(auto&& fn) const
    {
        Constexprifier{info_.isResistive(), info_.isHyperResistive()}(fn);
    }

    // block per tile, block threads stride each E component's box
    //  waits for the launch before returning: the per tile items live in this scope
    //   todo: async (no wait), the items must then outlive the launch, see .gpu_todo
    template<typename Field, typename VecField>
    void on_tiles_gpu([[maybe_unused]] Field const& n, [[maybe_unused]] VecField const& Ve,
                      [[maybe_unused]] Field const& Pe, [[maybe_unused]] VecField const& B,
                      [[maybe_unused]] VecField const& J, [[maybe_unused]] VecField& Enew) const
    {
#if PHARE_HAVE_GPU
        using Launcher = gpu::FieldTileLauncher</*sync=*/true>;

        auto const make_item = [&](std::size_t const ti) {
            auto const& layout = Enew[0][ti].layout();
            using Layout_t     = std::decay_t<decltype(layout)>;
            auto const nt      = tile_at(n, ti);
            auto const Pet     = tile_at(Pe, ti);
            auto const Vet     = tile_at(Ve, ti);
            auto const Bt      = tile_at(B, ti);
            auto const Jt      = tile_at(J, ti);
            auto Et            = tile_at(Enew, ti);
            auto const boxes = for_N_make_array<3>([&](auto c) { return box_for(layout, Et[c]); });
            auto const gradPe = grad(layout, Pet);
            using Box_t       = std::decay_t<decltype(boxes)>::value_type;
            using Item_t      = OhmTile<Ohm<Layout_t>, std::decay_t<decltype(gradPe)>,
                                        std::decay_t<decltype(nt)>, decltype(Et), Box_t>;
            return Item_t{Ohm<Layout_t>{info_, layout}, gradPe, nt, Vet, Bt, Jt, Et, boxes};
        };
        using Item_t = decltype(make_item(0));

        auto const n_tiles = tile_count(Enew);
        gpu::Vec_t<Item_t> items;
        items.reserve(n_tiles);
        for (std::size_t ti = 0; ti < n_tiles; ++ti)
            items.emplace_back(make_item(ti));
        auto* const ip = items.data();

        with_flags([&]<bool isResistive, bool isHyperResistive>() {
            Launcher{n_tiles}([=] _PHARE_ALL_FN_() mutable {
                auto& it        = ip[Launcher::tile_idx()];
                auto const eval = [&]<auto c>(auto const& box, auto const& gradPe_c) {
                    Launcher::for_each_point(box, [&](auto const& ijk) {
                        it.ohm.template operator()<c, isResistive, isHyperResistive>(
                            ijk, it.n, it.Ve, gradPe_c, it.B, it.J, it.E);
                    });
                };
                eval.template operator()<Component::X>(it.boxes[0], it.gradPe.x);
                eval.template operator()<Component::Y>(it.boxes[1], it.gradPe.y);
                eval.template operator()<Component::Z>(it.boxes[2], it.gradPe.z);
            });
        });
#else
        throw std::runtime_error("OhmTransformer: GPU_UNIFIED tiles need a GPU build");
#endif
    }

    OhmInfo info_;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_NUMERICS_OHM_OHM_TRANSFORMER_HPP */
