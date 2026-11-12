#ifndef PHARE_CORE_UTILITIES_KERNELS_REGION_CHUNKS_HPP
#define PHARE_CORE_UTILITIES_KERNELS_REGION_CHUNKS_HPP

// device launches over regions (a patch's box, or a tile's) collected across a patch or a level:
//  each (region, component) box is split evenly into chunks, a block per chunk, threads stride it,
//  one launch for all. Even split: ceil(size / target) chunks of ceil(size / n) points, so no
//  chunk is a mostly idle remainder. The target sizes the grid to the device (SM count)

#include "core/def.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <stdexcept>

namespace PHARE::core
{

struct RegionChunk
{
    std::uint32_t item;  // index in the items list
    std::uint32_t comp;  // component of the item, its boxes[comp]
    std::uint32_t start; // first point index in that box
    std::uint32_t count;
};

// points per chunk for a launch over total_points: ~8 blocks per SM, within [1024, 16384]
inline std::size_t region_chunk_target([[maybe_unused]] std::size_t const total_points)
{
#if PHARE_HAVE_GPU
    std::size_t constexpr static blocks_per_sm = 8, lo = 1 << 10, hi = 1 << 14;
    return std::clamp<std::size_t>(total_points / (blocks_per_sm * gpu::sm_count()), lo, hi);
#else
    return 1 << 12;
#endif
}

inline void append_region_chunks(std::uint32_t const item, std::uint32_t const comp,
                                 std::size_t const size, std::size_t const target, auto& chunks)
{
    if (size == 0)
        return;
    auto const n     = (size + target - 1) / target;
    auto const chunk = (size + n - 1) / n;
    for (std::size_t start = 0; start < size; start += chunk)
        chunks.push_back(RegionChunk{item, comp, static_cast<std::uint32_t>(start),
                                     static_cast<std::uint32_t>(std::min(chunk, size - start))});
}

// fn(item, comp, point) on the device for every point of every items[i].boxes[comp], one launch
//  items: a gpu::Vec_t of structs with std::array boxes, device copyable
template<typename Items, typename Fn>
void run_region_chunks([[maybe_unused]] Items& items, [[maybe_unused]] Fn const fn)
{
#if PHARE_HAVE_GPU
    if (items.empty())
        return;

    std::size_t total = 0;
    for (auto const& item : items)
        for (auto const& box : item.boxes)
            total += box.size();
    auto const target = region_chunk_target(total);

    static thread_local gpu::Vec_t<RegionChunk> chunks;
    chunks.clear();
    for (std::uint32_t i = 0; i < items.size(); ++i)
        for (std::uint32_t c = 0; c < items[i].boxes.size(); ++c)
            append_region_chunks(i, c, items[i].boxes[c].size(), target, chunks);

    using Launcher = gpu::FieldTileLauncher</*sync=*/true>;
    auto* const ip = items.data();
    auto* const cp = chunks.data();
    Launcher{chunks.size()}([=] _PHARE_ALL_FN_() mutable {
        auto const chunk = cp[Launcher::tile_idx()];
        auto& item       = ip[chunk.item];
        auto const& box  = item.boxes[chunk.comp];
        Launcher::for_each_index(chunk.count, [&](auto const i) {
            fn(item, chunk.comp, point_at(box, chunk.start + i));
        });
    });
#else
    throw std::runtime_error("run_region_chunks: needs a GPU build");
#endif
}

} // namespace PHARE::core

#endif /* PHARE_CORE_UTILITIES_KERNELS_REGION_CHUNKS_HPP */
