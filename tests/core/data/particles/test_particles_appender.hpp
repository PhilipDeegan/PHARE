#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP

#include <string>
#include <cassert>
#include <functional>

#include "phare_core.hpp"

#include "test_particles_converter.hpp"

namespace PHARE::core
{

// NOT FOR GPU
template<auto type, typename Src, typename Dst>
struct ParticlesAppender
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(Src::dimension == Dst::dimension);
    auto constexpr static src_alloc_mode = Src::alloc_mode;
    auto constexpr static dst_alloc_mode = Dst::alloc_mode;

    auto constexpr static src_layout_mode = Src::layout_mode;
    auto constexpr static dst_layout_mode = Dst::layout_mode;

    auto& operator()(Src const& src, Dst& dst)
    {
        for (auto const& p : src)
            dst.push_back(p);

        if constexpr (any_in(Dst::layout_mode, LayoutMode::AoSPC, LayoutMode::SoAPC))
            dst.template sync<2, type>();

        return dst;
    }
};



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
    // hax
    auto const ppc = [&]() { // assume constant for the moment
        for (auto const& bix : src.box())
            if (auto const n = src.nbr_particles_in(*bix); n > 0)
                return n;
        throw std::runtime_error("ppc 0");
    }();

    using Tmp
        = ParticleArray<Dst::dimension,
                        ParticleArrayInternals<Dst::dimension, LayoutMode::SoA, StorageMode::VECTOR,
                                               src_alloc_mode, Dst::impl>>;

    auto tmp = make_particles<Tmp>(layout);
    tmp.resize(src.size());
    ParticlePacker<Src>{src}.pack(tmp);
    assert(array_equals(src.delta(0), tmp.delta(0)));

    std::size_t tmp_start = 0;
    auto const tmp_tuple  = tmp.as_tuple();
    auto finish           = [&](auto& lix) {
        auto const size      = ppc;
        auto& dst_arr        = dst(lix);
        auto const curr_size = dst_arr.size();
        Dst::resize(dst_arr, curr_size + size);
        auto dst_tuple = dst_arr.as_tuple();

        for_N<std::tuple_size_v<decltype(dst_tuple)>>([&](auto vi) {
            auto const& tmp_vec = std::get<vi>(tmp_tuple);
            auto& vec           = std::get<vi>(dst_tuple);
            PHARE_WITH_MKN_GPU(
                mkn::gpu::copy(vec.data() + curr_size, tmp_vec.data() + tmp_start, size); //
            )
        });
        tmp_start += size;
    };

    auto const amrbox = dst.safe_box();
    auto const lclbox = dst.local_box();
    auto amr_it       = amrbox.begin();
    auto lcl_it       = lclbox.begin();

    for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
    {
        if (!array_equals(**amr_it, tmp.iCell(tmp_start)))
            continue;
        finish(*lcl_it);
    }

    dst.template sync<2, type>();
    assert(dst.size() == src.size());

    bool valid = false;
    amr_it     = amrbox.begin();
    lcl_it     = lclbox.begin();

    for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
    {
        if (!array_equals(**amr_it, tmp.iCell(0)))
            continue;
        valid = array_equals(src.delta(0), dst(*lcl_it).delta(0));
    }

    auto const& last = [&]() -> auto& {
        static_assert(std::is_same_v<decltype(type), ParticleType>);
        if constexpr (type == ParticleType::Domain)
            return dst(dst.local_cell(dst.box().upper));
        else if (type == ParticleType::Ghost)
            return dst(dst.local_cell(dst.ghost_box().upper));
    }();
    assert(last.size() == ppc);
    assert(valid);
}




template<>
template<auto type, typename Src, typename Dst, typename GridLayout>
void PAppender<LayoutMode::AoSMapped, AllocatorMode::CPU, LayoutMode::AoSPC,
               AllocatorMode::GPU_UNIFIED>::operator()(Src const& src, Dst& dst,
                                                       GridLayout const& layout)
{
    // hax
    auto const ppc = [&]() { // assume constant for the moment
        for (auto const& bix : src.box())
            if (auto const n = src.nbr_particles_in(*bix); n > 0)
                return n;
        throw std::runtime_error("ppc 0");
    }();

    std::size_t src_start = 0;
    auto finish           = [&](auto& lix) {
        auto const size      = ppc;
        auto& dst_arr        = dst(lix);
        auto const curr_size = dst_arr.size();
        Dst::resize(dst_arr, curr_size + size);
        PHARE_WITH_MKN_GPU(
            mkn::gpu::copy(dst_arr.data() + curr_size, src.data() + src_start, size); //
        )
        src_start += size;
    };

    auto const amrbox = dst.safe_box();
    auto const lclbox = dst.local_box();
    auto amr_it       = amrbox.begin();
    auto lcl_it       = lclbox.begin();

    for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
    {
        if (!array_equals(**amr_it, src.iCell(src_start)))
            continue;
        finish(*lcl_it);
    }

    dst.template sync<2, type>();
    assert(dst.size() == src.size());

    bool valid = false;
    amr_it     = amrbox.begin();
    lcl_it     = lclbox.begin();

    for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
    {
        if (!array_equals(**amr_it, src.iCell(0)))
            continue;
        valid = array_equals(src.delta(0), dst(*lcl_it).delta(0));
    }

    auto const& last = [&]() -> auto& {
        static_assert(std::is_same_v<decltype(type), ParticleType>);
        if constexpr (type == ParticleType::Domain)
            return dst(dst.local_cell(dst.box().upper));
        else if (type == ParticleType::Ghost)
            return dst(dst.local_cell(dst.ghost_box().upper));
    }();
    assert(last.size() == ppc);
    assert(valid);
}



template<auto type, typename Src, typename Dst, typename GridLayout>
void append(Src const& src, Dst& dst, GridLayout const& layout)
{
    PAppender<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>{}
        .template operator()<type>(src, dst, layout);
}



} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_APPENDER_HPP */
