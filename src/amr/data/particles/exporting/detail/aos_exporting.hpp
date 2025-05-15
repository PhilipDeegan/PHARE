
#ifndef PHARE_AMR_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER
#define PHARE_AMR_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER


#include "core/def.hpp"

// #include "core/utilities/box/box.hpp"
// #include "core/data/particles/particle_array_def.hpp"
// #include "core/data/particles/particle_array_appender.hpp"
// #include "core/data/particles/particle_array_partitioner.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/exporting/detail/def_exporting.hpp"


#include <stdexcept>


namespace PHARE::amr
{

using enum core::LayoutMode;
using enum AllocatorMode;


template<>
template<typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesExporter<AoS, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box, Fn0 fn0,
                                             Fn1 fn1)
{
    std::uint16_t constexpr static N = 256;

    auto const coarseDstBox = coarsen_box(box);

    using ArrayParticleArray = typename Src::template array_type<N>;
    using SpanParticleArray  = Src::Span_t;

    std::uint16_t big_buffer_cnt = 0;
    ArrayParticleArray big_buffer;
    SpanParticleArray span{big_buffer};

    auto const send = [&]() {
        assert(big_buffer_cnt <= N);
        span.resize(big_buffer_cnt);
        append_particles<core::ParticleType::Domain>(span, dst);
        big_buffer_cnt = 0;
    };

    for (auto const& particle : src)
    {
        if (not isIn(particle, coarseDstBox))
            continue;

        if (not isIn(fn0(particle), box))
            continue;

        auto&& [p_count, buffer] = fn1(particle, box);
        if (big_buffer_cnt + p_count > N)
            send();

        std::copy(buffer.data(), buffer.data() + p_count, big_buffer.data() + big_buffer_cnt);
        big_buffer_cnt += p_count;
    }

    if (big_buffer_cnt)
        send();
}

template<>
template<typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesExporter<AoSMapped, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                   Fn0 fn0, Fn1 fn1)
{
    // auto const old_size = dst.size();
    ParticlesExporter<AoS, CPU>{}(src, dst, box, fn0, fn1); // we don't use cellmap here
    // dst.map_particles(old_size);
}


template<>
template<typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesExporter<AoSTS, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box, Fn0 fn0,
                                               Fn1 fn1)
{
    auto const& dstBox = dst.box();

    auto const fineDstOverlapOpt = box * dstBox;
    if (!fineDstOverlapOpt)
        return;

    auto const& fineDstOverlap = *fineDstOverlapOpt;
    auto const coarseDstBox    = coarsen_box(box);

    for (auto const& src_tile : src()) // !expensive! covers ghost box
    {
        auto const overlap_opt = coarseDstBox * src_tile;
        if (not overlap_opt)
            continue;

        auto const overlap        = *overlap_opt;
        auto const fineOverlapBox = refine_box(overlap);

        auto const src_dst_overlap_opt = fineOverlapBox * fineDstOverlap;
        if (!src_dst_overlap_opt)
            continue;

        auto const& src_dst_overlap = *src_dst_overlap_opt;
        auto& dst_tile              = *dst().at(src_dst_overlap.lower);
        auto const dst_tile_overlap = src_dst_overlap * dst_tile;
        assert(dst_tile_overlap);
        ParticlesExporter<AoS, CPU>{}(src_tile(), dst_tile(), *dst_tile_overlap, fn0, fn1);


        // maybe needed?
        // for (auto& link_opt : dst_tile.links())
        //     if (link_opt and)
        //         ParticlesExporter<AoS, CPU>{}(tile(), (*link_opt)(), box, fn0, fn1);
        // PHARE_LOG_LINE_SS(dst_tile().size());
    }

    PHARE_DEBUG_DO({
        for (auto const& tile : dst())
            for (auto const& p : tile())
                if (!isIn(p, tile))
                    throw std::runtime_error("no");
    })
}


template<> // slow
template<typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesExporter<AoS, GPU_UNIFIED>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                     Fn0 fn0, Fn1 fn1)
{
    ParticlesExporter<AoS, CPU>{}(src, dst, box, fn0, fn1);
}


template<> // slow
template<typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                       Fn0 fn0, Fn1 fn1)
{
    ParticlesExporter<AoSTS, CPU>{}(src, dst, box, fn0, fn1);
}




} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER */
