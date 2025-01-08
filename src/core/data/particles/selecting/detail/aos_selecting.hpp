
#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING


#include "core/utilities/memory.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


// AoS, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                         box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoS, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                         box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

// AoS, GPU_UNIFIED

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                 box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoS, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                 box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

// AoSPC, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                           box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto const& p : sv)
            dv.emplace_back(p);
    }
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoSPC, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                           box_t const& box, Shift&& fn)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - fn); // BEWARE
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto p : sv)
        {
            p.iCell() = (Point{p.iCell()} + fn).toArray();
            dv.emplace_back(p);
        }
    }
}

//

// AoSPC, GPU_UNIFIED

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                   box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto const& p : sv)
            dv.emplace_back(p);
    }
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                   box_t const& box, Shift&& shift)
{ // unoptimized
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - shift);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto const& pit : sv)
        {
            auto p    = pit.copy();
            p.iCell() = (Point{p.iCell()} + shift).toArray();
            dv.emplace_back(p);
        }
    }
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING */
