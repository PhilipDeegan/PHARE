#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING

#include "core/def/phare_config.hpp"
#include "core/utilities/memory.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"
// #include "core/data/particles/particle_array_converter.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"


namespace PHARE::core
{

using enum LayoutMode;

// CPU-AoSMapped to GPU_UNIFIED-SoAPC
template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void ParticlesAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoAPC,
                       AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                               GridLayout const& layout)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    using Tmp // make chunked
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
void ParticlesAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoATS,
                       AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                               GridLayout const& /*layout*/)
{
    static constexpr std::uint8_t N = 128;

    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<dst_alloc_mode>(a.data() + particles.size(), b.data(), S);
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
void ParticlesAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoATS,
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
void ParticlesAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::SoAVXTS,
                       AllocatorMode::CPU>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& /*layout*/)
{
    static constexpr std::uint8_t N = 128;

    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<dst_alloc_mode>(a.data() + particles.size(), b.data(), S);
        });
        particles.resize(particles.size() + S);
    };

    std::int32_t src_start = 0;

    auto do_copy = [&](auto&&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        for (std::size_t pidx = 0; pidx < S; ++pidx)
            tmp.assign(src[pidx + src_start], pidx);
        copy(particles, tmp, i, S);
        src_start += S;
    };

    SoAVXArray<GridLayout::dimension, N> tmp;

    auto finish = [&](auto const lix, auto const size) {
        auto& particles       = dst(lix);
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

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING */
