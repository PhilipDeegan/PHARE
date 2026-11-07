// IWYU pragma: private, include "core/data/particles/appending/particles_appending.hpp"
#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING


#include "core/utilities/memory.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"


#include <iterator>

namespace PHARE::core
{

using LM = LayoutMode;
using AM = AllocatorMode;


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    // reserve_flat's upfront estimate is average-case, not a hard bound - reserving here on
    // every chunk would force an exact-size reallocation instead of amortized push_back
    // growth for any shortfall.
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoSMapped, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}




template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}

template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSCMTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSCMTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPCTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSPCTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}




template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoSPC, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPC, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
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

    for_each_point(*overlap, [&](auto const& bix) {
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);
    });

    dst.template on_appended<type>();
}


// into one AoSPCTS tile (AoSPC), e.g. the refined particle chunks of the AoSPCTS GPU refiner
//  (forwards to the CPU one, src is a span over a host buffer)
template<> // slow - unified memory, CPU path
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::GPU_UNIFIED, LM::AoSPC, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoS, GPU_UNIFIED, AoSPC, GPU_UNIFIED>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}


template<> // slow - unified memory, CPU path
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPCTS, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPCTS, AM::CPU>{}.template operator()<type>(
        src, dst);
}


// src is host memory: staged contiguously into managed memory, then counted and scattered into
//  the dst tiles on the device, see TileSetVector::append_on_device. Writing dst tiles from the
//  host would fault their (device resident) pages back
template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSTS, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSTS, GPU_UNIFIED>::operator()");

    if (src.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    using Particle_t = typename Dst::Particle_t;
    gpu::Vec_t<Particle_t> staged(src.size());
    std::copy(src.data(), src.data() + src.size(), staged.data());
    dst.template append_on_device<type>(staged.data(), staged.size());
#else
    for (auto const& p : src)
        dst.emplace_back(p);
#endif
}



} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING */
