
#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/utilities/memory.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/converting/detail/def_converting.hpp"

namespace PHARE::core
{



template<typename Dst, typename... Args>
void make_return_object(Args&&... args)
{
    using enum StorageMode;

    static_assert(any_in(Dst::storage_mode, VECTOR, ARRAY));
    if constexpr (Dst::storage_mode == StorageMode::VECTOR)
        return make_particles(args...);
    else if constexpr (Dst::storage_mode == StorageMode::ARRAY)
        return Dst{};
};

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

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING */