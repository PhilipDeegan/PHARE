#ifndef PHARE_CORE_DATA_FIELD_FIELD_OVERLAPS_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_OVERLAPS_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/field/field_box.hpp"
#include "core/data/grid/grid_tiles.hpp"


#include <cstdint>
#include <vector>
#include <cstddef>
#include <cstring>



namespace PHARE::core
{


template<auto opts>
struct FieldTileOverlaps
{
    // when you're sure it should exist
    static auto& getInnerTileOverlaps(int const lvl, auto const& layout)
    {
        return levels.at(lvl).at(to_string(layout.AMRBox())).tiles;
    }

    // when you're not sure it should exist
    static auto& getOrCreateInnerTileOverlaps(int const lvl, auto const& layout, auto const& field);

    static auto build_inner_tile_overlaps(auto const& field);
    static void sync_inner_ghosts(auto& field, auto const& overlaps);
    static void reset(int lvl) { levels[lvl].clear(); }
    static void reset() { levels.clear(); }


    struct Overlap
    {
        Overlap(auto _src, auto const d, auto const s)
            : src{_src}
            , lcl_dst{d}
            , lcl_src{s}
        {
        }

        basic::Field<opts> const* const src;
        Box<std::uint32_t, opts.dimension> lcl_dst, lcl_src;
    };

    struct Level
    {
        struct Patch
        {
            struct Tile
            {
                std::vector<Overlap> overlaps{};
            };

            std::vector<Tile> tiles{};
        };

        std::map<std::string, Patch> patches{};
    };

    static inline std::map<int, Level> levels{};
};

template<auto opts>
auto& FieldTileOverlaps<opts>::getOrCreateInnerTileOverlaps( //
    int const lvl, auto const& layout, auto const& field)
{
    if (!levels.count(lvl))
        levels.try_emplace(lvl);

    auto& level           = levels.at(lvl);
    auto const& patch_key = to_string(layout.AMRBox());
    if (!level.patches.count(patch_key))
        level.patches.try_emplace(patch_key);

    auto& patch = level.patches.at(patch_key);
    if (!patch.tiles.size())
        patch.tiles = build_inner_tile_overlaps(field);

    return patch.tiles;
}

template<auto opts>
auto FieldTileOverlaps<opts>::build_inner_tile_overlaps(auto const& field)
{
    using TileOverlap = FieldTileOverlaps<opts>::Level::Patch::Tile;

    std::vector<TileOverlap> tiles(field().size());

    for (std::size_t ti0 = 0; ti0 < field().size() - 1; ++ti0)
        for (std::size_t ti1 = ti0 + 1; ti1 < field().size(); ++ti1)
        {
            auto& t0 = field()[ti0];
            auto& t1 = field()[ti1];

            if (auto const t0_overlap = t0.ghost_box() * t1.field_box())
                tiles[ti0].overlaps.emplace_back(
                    &t1(), as_unsigned(shift(*t0_overlap, t0.ghost_box().lower * -1)),
                    as_unsigned(shift(*t0_overlap, t1.ghost_box().lower * -1)));
            else
                continue; // second can't overlap if first doesn't

            if (auto const t1_overlap = t1.ghost_box() * t0.field_box())
                tiles[ti1].overlaps.emplace_back(
                    &t0(), as_unsigned(shift(*t1_overlap, t1.ghost_box().lower * -1)),
                    as_unsigned(shift(*t1_overlap, t0.ghost_box().lower * -1)));
        }


    return tiles;
}

template<auto opts>
void FieldTileOverlaps<opts>::sync_inner_ghosts(auto& field, auto const& overlaps_per_tile)
{
    using value_type = decltype(opts)::value_type;

    for (std::size_t i = 0; i < field().size(); ++i)
        for (auto const& overlap : overlaps_per_tile[i].overlaps)
            operate_on_fields<Equals<value_type>>(
                FieldBox{field()[i], field()[i].layout(), overlap.lcl_dst},
                FieldBox{*overlap.src, field()[i].layout(), overlap.lcl_src});
}


} // namespace PHARE::core


#endif
