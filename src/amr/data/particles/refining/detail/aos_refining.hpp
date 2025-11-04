
#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER


#include "core/def.hpp"

// #include "core/utilities/box/box.hpp"
// #include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
// #include "core/data/particles/particle_array_partitioner.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/refining/detail/def_refining.hpp"




namespace PHARE::amr
{

using enum core::LayoutMode;
using enum AllocatorMode;


template<>
template<auto type, typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesRefiner<AoS, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box, Fn0 fn0,
                                            Fn1 fn1)
{
    // PHARE_LOG_LINE_SS(box);
    std::uint16_t constexpr static N = 256;

    auto const splitBox     = grow(box, 2); // Splitter::maxCellDistanceFromSplit() ?
    auto const coarseDstBox = coarsen_box(splitBox);

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

        if (not isIn(fn0(particle), splitBox))
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
template<auto type, typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesRefiner<AoSMapped, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                  Fn0 fn0, Fn1 fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<type>(src, dst, box, fn0, fn1);
}


template<>
template<auto type, typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesRefiner<AoSTS, CPU>::operator()(Src const& src, Dst& dst, Box_t const& box, Fn0 fn0,
                                              Fn1 fn1)
{
    // PHARE_LOG_LINE_SS(box);
    auto const& dstBox  = dst.box();
    auto const splitBox = grow(box, 1); // Splitter::maxCellDistanceFromSplit() ?

    auto const ghost_nbr = (dst.ghost_box().shape()[0] - dst.box().shape()[0]) / 2;

    auto const coarseDstBox = coarsen_box(splitBox);

    for (auto const& src_tile : src())
    {
        if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
            continue;

        if constexpr (type == core::ParticleType::Domain)
        {
            for (auto& dst_tile : dst())
                if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                    if (auto const overlap = box * dst_tile)
                        ParticlesRefiner<AoS, CPU>{}.template operator()<type>(
                            src_tile(), dst_tile(), *overlap, fn0, fn1);
        }
        else if constexpr (type == core::ParticleType::Ghost)
        {
            for (auto& dst_tile : dst()) // only on dst patch borders
                if (auto const growbox = grow(dst_tile, ghost_nbr);
                    *(dst.box() * growbox) != growbox)
                    if (auto const overlap = box * growbox)
                        ParticlesRefiner<AoS, CPU>{}.template operator()<type>(
                            src_tile(), dst_tile(), *overlap, fn0, fn1);
        }
        else
        {
            assert(false);
        }
    }
}



template<> // slow
template<auto type, typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesRefiner<AoS, GPU_UNIFIED>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                    Fn0 fn0, Fn1 fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<type>(src, dst, box, fn0, fn1);
}


template<> // slow
template<auto type, typename Src, typename Dst, typename Box_t, typename Fn0, typename Fn1>
void ParticlesRefiner<AoSTS, GPU_UNIFIED>::operator()(Src const& src, Dst& dst, Box_t const& box,
                                                      Fn0 fn0, Fn1 fn1)
{
    ParticlesRefiner<AoSTS, CPU>{}.template operator()<type>(src, dst, box, fn0, fn1);
}




} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER */
