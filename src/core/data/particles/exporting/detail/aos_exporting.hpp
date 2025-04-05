
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER


#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"
#include <optional>


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


template<>
template<typename Src, typename Dst, typename Box_t, typename Fn0>
void ParticlesExporter<AoSMapped, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                   Fn0 fn0)
{
    constexpr static std::uint16_t N = 256;

    using ArrayParticleArray = typename Src::template array_type<N>;
    using SpanParticleArray  = Src::Span_t;

    std::uint16_t big_buffer_cnt = 0;
    ArrayParticleArray big_buffer;
    SpanParticleArray span{big_buffer};

    auto const send = [&]() {
        assert(big_buffer_cnt <= N);
        span.resize(big_buffer_cnt);
        append_particles<ParticleType::Domain>(span, dst);
        big_buffer_cnt = 0;
    };

    for (auto const& particle : src)
    {
        if (not isIn(particle, box))
            continue;

        auto&& [p_count, buffer] = fn0(particle);
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
void ParticlesExporter<AoS, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box, Fn0 fn0,
                                             Fn1 fn1)
{
    std::uint16_t constexpr static N   = 256;
    std::size_t constexpr static ratio = 2; // to do - get somewhere else

    auto const coarseDstBox = box / ratio;

    using ArrayParticleArray = typename Src::template array_type<N>;
    using SpanParticleArray  = Src::Span_t;

    std::uint16_t big_buffer_cnt = 0;
    ArrayParticleArray big_buffer;
    SpanParticleArray span{big_buffer};

    auto const send = [&]() {
        assert(big_buffer_cnt <= N);
        span.resize(big_buffer_cnt);
        append_particles<ParticleType::Domain>(span, dst);
        big_buffer_cnt = 0;
    };

    for (auto const& particle : src)
    {
        if (not isIn(particle, coarseDstBox))
            continue;

        if (not isIn(fn0(particle), box))
            continue;

        auto&& [p_count, buffer] = fn1(particle);
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
    std::size_t constexpr ratio = 2; // to do - get somewhere else

    auto const& dstBox           = dst.ghost_box();
    auto const fineDstOverlapOpt = box * dstBox;
    if (!fineDstOverlapOpt)
        return;

    auto const& fineDstOverlap = *fineDstOverlapOpt;
    auto const coarseDstBox    = box / ratio;

    for (auto const& src_tile : src()) // !expensive!
    {
        auto const overlap_opt = coarseDstBox * src_tile;
        if (not overlap_opt)
            continue;

        auto const overlap        = *overlap_opt;
        auto const fineOverlapBox = overlap * ratio;

        auto const src_dst_overlap_opt = fineOverlapBox * fineDstOverlap;
        if (!src_dst_overlap_opt)
            continue;

        auto const& src_dst_overlap = *src_dst_overlap_opt;
        auto& dst_tile              = *dst().at(src_dst_overlap.lower);

        ParticlesExporter<AoS, CPU>{}(src_tile(), dst_tile(), box, fn0, fn1);


        // maybe needed?
        // for (auto& link_opt : dst_tile.links())
        //     if (link_opt and)
        //         ParticlesExporter<AoS, CPU>{}(tile(), (*link_opt)(), box, fn0, fn1);
        // PHARE_LOG_LINE_SS(dst_tile().size());
    }
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


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER */
