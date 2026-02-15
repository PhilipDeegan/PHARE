
#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING


// #include "core/utilities/memory.hpp"
#include "core/logger.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"


#include <stdexcept>

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


//

// AoS, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

// AoSTS, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    auto on_overlap = [](auto& a, auto& b, auto f) {
        if (auto const overlap = a * b)
            f(*overlap);
    };

    auto const per_tile = [&](auto& src_tile, auto& dst_tile) { //
        auto& dst_particles     = dst_tile();
        auto& src_particles     = src_tile();
        auto const old_dst_size = dst_particles.size();

        on_overlap(box, *src_tile, [&](auto const& o0) {
            on_overlap(o0, *dst_tile, [&](auto const& overlap) {
                if (overlap == src_tile) // complete
                {
                    dst_particles.resize(old_dst_size + src_particles.size());
                    mem::copy<CPU>(dst_particles.data() + old_dst_size, src_particles.data(),
                                   src_particles.size());
                }
                else
                {
                    auto const in_dst = partition_particles(src_particles, overlap).size();
                    dst_particles.resize(old_dst_size + in_dst);
                    mem::copy<CPU>(dst_particles.data() + old_dst_size, src_particles.data(),
                                   in_dst);
                }
            });
        });

    };

    auto& mutable_src = const_cast<SrcParticles&>(src); // :(
    traverse_tilesets_overlap(mutable_src /*.views()*/, dst /*()*/, box, per_tile);
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoSTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    using enum LayoutMode;
    static_assert(any_in(DstParticles::layout_mode, AoS, AoSTS), "Otherwise not implemented!");

    auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
    auto const copy        = [&](auto& s, auto& d, auto const& b) {
        auto const in_dst = partition_particles(s, b);
        if (in_dst.size() == 0)
            return;

        auto const old_dst_size = d.size();
        d.resize(d.size() + in_dst.size());

        mem::copy<CPU>(d.data() + old_dst_size, s.data(), in_dst.size());
        if (any(shifter, [](auto const e) { return e != 0; }))
            for (std::size_t i = old_dst_size; i < d.size(); i++)
                d[i].iCell() = *((shifter) + d[i].iCell());
    };

    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        if (auto const src_overlap = grow(src_tile, ghost_cells) * box)
        {
            if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
            {
                for (auto& dst_tile : dst())
                {
                    auto const dst_gbox = shift(dst_tile, shifter * -1);
                    if (auto const overlap = dst_gbox * *src_overlap)
                        copy(src_tile(), dst_tile(), *overlap);
                }
            }
            else
            {
                copy(src_tile(), dst, *src_overlap);
            }
        }
    }

    if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
        dst.template sync<2>();
}


//

// AoS, GPU_UNIFIED

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

// AoSPC, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
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
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
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
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
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
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shift)
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

//

// AoSMapped, CPU

template<>
template<typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    src.export_particles(box, dst);
}

template<>
template<typename SrcParticles, typename DstParticles, typename box_t, typename Shift>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    auto const offseter = [&](auto const& particle) {
        auto shiftedParticle{particle};
        for (std::size_t idir = 0; idir < SrcParticles::dimension; ++idir)
            shiftedParticle.iCell()[idir] += shifter[idir];
        return shiftedParticle;
    };

    src.export_particles(box, dst, offseter);
}

//


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING */
