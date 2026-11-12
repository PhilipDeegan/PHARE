#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOAVX_CHUNKED_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOAVX_CHUNKED_APPENDING

// SoAVXChunked from / to AoS (pack and refine staging) and itself, one particle at a time,
//  GPU_UNIFIED on the host. The destination's offsets are rebuilt by the caller's on_appended

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"

namespace PHARE::core
{

template<typename Src, typename Dst>
void append_soavx_chunked(Src const& src, Dst& dst, std::size_t const start,
                          std::size_t const end)
{
    // refine streams in small chunks on top of reserve_flat's average case estimate: an exact
    //  reserve per chunk reallocates and copies every column each time (quadratic, and on
    //  GPU_UNIFIED the umpire pool grew to >24GB), so only grow when short, with 5% headroom
    double constexpr static over_alloc = 1.05;
    if (auto const needed = dst.size() + end - start; needed > dst.capacity())
        dst.reserve(static_cast<std::size_t>(needed * over_alloc));
    for (std::size_t i = start; i < end; ++i)
        if constexpr (is_soavx(Src::layout_mode)) // columns: a copy, no particle reference
            dst.emplace_back(src.copy(i));
        else
            dst.emplace_back(src[i]);
}

template<auto src_alloc_mde, auto dst_alloc_mde>
struct ParticlesAppender<LayoutMode::AoS, src_alloc_mde, LayoutMode::SoAVXChunked, dst_alloc_mde>
{
    template<auto type, typename Src, typename Dst>
    void operator()(Src const& src, Dst& dst)
    {
        append_soavx_chunked(src, dst, start, end);
    }

    std::size_t const start;
    std::size_t const end;
};

template<auto src_alloc_mde, auto dst_alloc_mde> // AoSMapped is AoS + CellMap: same per particle
struct ParticlesAppender<LayoutMode::AoSMapped, src_alloc_mde, LayoutMode::SoAVXChunked,
                         dst_alloc_mde>
{
    template<auto type, typename Src, typename Dst>
    void operator()(Src const& src, Dst& dst)
    {
        append_soavx_chunked(src, dst, start, end);
    }

    std::size_t const start;
    std::size_t const end;
};

template<auto src_alloc_mde, auto dst_alloc_mde>
struct ParticlesAppender<LayoutMode::SoAVXChunked, src_alloc_mde, LayoutMode::SoAVXChunked,
                         dst_alloc_mde>
{
    template<auto type, typename Src, typename Dst>
    void operator()(Src const& src, Dst& dst)
    {
        append_soavx_chunked(src, dst, start, end);
    }

    std::size_t const start;
    std::size_t const end;
};

template<auto src_alloc_mde, auto dst_alloc_mde>
struct ParticlesAppender<LayoutMode::SoAVXChunked, src_alloc_mde, LayoutMode::AoS, dst_alloc_mde>
{
    template<auto type, typename Src, typename Dst>
    void operator()(Src const& src, Dst& dst)
    {
        append_soavx_chunked(src, dst, start, end);
    }

    std::size_t const start;
    std::size_t const end;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOAVX_CHUNKED_APPENDING */
