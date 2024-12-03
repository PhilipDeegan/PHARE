#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP


#include "core/def.hpp"
#include "core/utilities/memory.hpp"
#include "core/data/particles/particle_packer.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cassert>

namespace PHARE::core
{

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct PAppender
{
    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<auto type, typename Src, typename Dst, typename GridLayout>
    void operator()(Src const& src, Dst& dst, GridLayout const& layout);
};


// CPU-AoSMapped to GPU_UNIFIED-SoAPC
template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoAPC,
               AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& layout)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    using Tmp
        = ParticleArray<Dst::dimension,
                        ParticleArrayInternals<Dst::dimension, LayoutMode::SoA, StorageMode::VECTOR,
                                               src_alloc_mode, Dst::impl>>;

    auto tmp = make_particles<Tmp>(layout);
    tmp.resize(src.size());
    ParticlePacker<Src>{src}.pack(tmp);

    std::size_t tmp_start = 0;
    auto const tmp_tuple  = tmp.as_tuple();
    auto finish           = [&](auto const lix, auto const size) {
        auto& dst_arr        = dst(lix);
        auto const curr_size = dst_arr.size();
        Dst::resize(dst_arr, curr_size + size);
        auto dst_tuple = dst_arr.as_tuple();
        for_N<std::tuple_size_v<decltype(dst_tuple)>>([&](auto vi) {
            auto const& tmp_vec = std::get<vi>(tmp_tuple);
            auto& vec           = std::get<vi>(dst_tuple);
            gpu::copy(vec.data() + curr_size, tmp_vec.data() + tmp_start, size);
        });
        tmp_start += size;
    };

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}




template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::AoSPC,
               AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        auto& dst_arr        = dst(lix);
        auto const curr_size = dst_arr.size();
        Dst::resize(dst_arr, curr_size + size);
        gpu::copy(dst_arr.data() + curr_size, src.data() + src_start, size);
        src_start += size;
    };

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}


template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::AoSTS,
               AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        auto& particles = dst(lix);
        mem::copy<AllocatorMode::GPU_UNIFIED>(particles.data() + particles.size(),
                                                         src.data() + src_start, size);
        particles.resize(particles.size() + size);
        src_start += size;
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    dst.reset_views();

    for (auto const bix : *overlap)
        if (auto const size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}

template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::AoSTS,
               AllocatorMode::CPU>::operator()(Src const& src, Dst& dst,
                                               GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        auto& dst_arr = dst(lix);
        dst_arr.append(src, src_start, size);
        src_start += size;
    };

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}




template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoATS,
               AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    static constexpr std::uint8_t N = 128;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<AllocatorMode::GPU_UNIFIED>(a.data() + particles.size(), b.data(), S);
        });
        particles.resize(particles.size() + S);
    };

    std::int32_t src_start = 0;
    auto do_copy           = [&](auto&&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        for (std::size_t pidx = 0; pidx < S; ++pidx)
            tmp.assign(src[pidx + src_start], pidx);
        copy(particles, tmp, i, S);
        src_start += S;
    };

    SoAArray<GridLayout::dimension, N> tmp;

    auto finish = [&](auto const lix, auto const size) {
        assert(dst().at(lix));
        auto& tile            = *dst().at(lix);
        auto& particles       = tile();
        auto const full_count = size / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
            do_copy(particles, tmp, i, N);
        if (auto const remaining = size % N)
            do_copy(particles, tmp, i, remaining);
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    dst.reset_views();

    for (auto const bix : *overlap)
        if (auto const size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}

template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoATS,
               AllocatorMode::CPU>::operator()(Src const& src, Dst& dst,
                                               GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        dst(lix).append(src, src_start, size);
        src_start += size;
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}


template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoAVXTS,
               AllocatorMode::CPU>::operator()(Src const& src, Dst& dst,
                                               GridLayout const& /*layout*/)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        auto const end = src_start + size;
        auto& arr      = dst(lix); //.append(src, src_start, size);
        for (std::size_t i = src_start; i < end; ++i)
        {
            //
        }
        src_start += size;
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    Src const& src, Dst& dst, GridLayout const& /*layout*/)
{
    using enum LayoutMode;
    static_assert(src_layout_mde == dst_layout_mde);
    static_assert(src_alloc_mde == dst_alloc_mde);
    static_assert(any_in(src_layout_mde, AoS, SoA));

    auto const oldSize = dst.size();
    dst.resize(oldSize + src.size());

    if constexpr (src_layout_mde == AoS)
    {
        mem::copy<src_alloc_mde>(dst.data() + oldSize, src.data(), src.size());
    }
    else if constexpr (src_layout_mde == SoA)
    {
        auto dst_tuple = dst.as_tuple();
        auto src_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(dst_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<src_alloc_mde>(a.data() + oldSize, b, src.size());
        });
    }
    else
        throw std::runtime_error("no");

    // dst.resize(dst.size() + src.size());

    // ParticleArrayService::sync<2, type>(dst);
}


template<auto type, typename Src, typename Dst, typename GridLayout>
void append_particles(Src const& src, Dst& dst, GridLayout const& layout)
{
    PHARE_DEBUG_DO(int const old_size = dst.size();)

    std::string_view constexpr static FN_ID = "append_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    PAppender<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>{}
        .template operator()<type>(src, dst, layout);

    PHARE_DEBUG_DO(assert(dst.size() == old_size + src.size());)
}



} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP */
