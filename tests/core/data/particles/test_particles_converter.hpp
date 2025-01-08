#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP


#include "core/utilities/memory.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include <tuple>


namespace PHARE::core
{
using enum LayoutMode;

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesConverter
{
    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<auto type, typename Dst, typename Src, typename GridLayout>
    Dst operator()(Src const& src, GridLayout const& layout);

    std::size_t const start;
    std::size_t const end;
};



template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, AllocatorMode::CPU, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
        std::copy(tile().begin(), tile().end(), std::back_inserter(out));

    return out;
}


template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, AllocatorMode::GPU_UNIFIED, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.resize(src.size());
    std::size_t offset = 0;
    for (auto const& tile : src())
    {
        auto& particles = src(src.local_cell(tile.lower));
        gpu::copy(out.data() + offset, particles.data(), particles.size());
        offset += particles.size();
    }
    return out;
}



template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoATS, AllocatorMode::CPU, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    static constexpr std::uint8_t N = 128;

    auto copy = [](auto&... args) {
        auto const& [tmp, particles, i, S] = std::forward_as_tuple(args...);
        auto tmp_tuple                     = tmp.as_tuple();
        auto src_tuple                     = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(tmp_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<AllocatorMode::CPU>(a.data(), b.data() + N * i, S);
        });
    };

    SoAArray<GridLayout::dimension, N> tmp;
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
    {
        auto const& particles = tile();
        auto const full_count = particles.size() / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
        {
            copy(tmp, particles, i, N);
            out.append(tmp, 0, N);
        }
        if (auto const remaining = particles.size() % N)
        {
            copy(tmp, particles, i, remaining);
            out.append(tmp, 0, remaining);
        }
    }
    assert(src.size() == out.size());
    return out;
}


template<>
template<auto type, typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoATS, AllocatorMode::GPU_UNIFIED, AoS, AllocatorMode::CPU>::operator()(
    Src const& src, GridLayout const& layout)
{
    static constexpr std::uint8_t N = 128;

    auto copy = [](auto&... args) {
        auto const& [tmp, particles, i, S] = std::forward_as_tuple(args...);
        auto tmp_tuple                     = tmp.as_tuple();
        auto src_tuple                     = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(tmp_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<AllocatorMode::GPU_UNIFIED>(a.data(), b.data() + N * i, S);
        });
    };

    SoAArray<GridLayout::dimension, N> tmp;
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
    {
        auto const& particles = tile();
        auto const full_count = particles.size() / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
        {
            copy(tmp, particles, i, N);
            out.append(tmp, 0, N);
        }
        if (auto const remaining = particles.size() % N)
        {
            copy(tmp, particles, i, remaining);
            out.append(tmp, 0, remaining);
        }
    }
    assert(src.size() == out.size());
    return out;
}



template<auto type, typename Dst, typename Src, typename... Args>
auto static convert_particles(Src const& src, Args&&... args)
{
    using Converter
        = ParticlesConverter<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>;

    std::string_view constexpr static FN_ID = "convert_particles,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    auto const& [layout, start, end] = std::forward_as_tuple(args...);
    return Converter{start, end}.template operator()<type, Dst>(src, layout);
}

template<auto type, typename Dst, typename Src, typename GridLayout>
auto static convert_particles(Src const& src, GridLayout const& layout)
{
    convert_particles<type, Dst>(src, layout, 0, src.size());
}

template<auto type, typename Dst, typename Src, typename GridLayout>
auto static convert_particles_and_sort(Src const& src, GridLayout const& layout)
{
    std::string_view constexpr static FN_ID = "convert_particles_and_sort,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    auto out = convert_particles<type, Dst>(src, layout);
    sort(out, layout.AMRBox());
    return out;
}

} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP */
