#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING

#include "core/utilities/memory.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"
#include "core/data/particles/appending/detail/aos_appending.hpp"
#include "core/data/particles/appending/detail/soa_appending.hpp"

namespace PHARE::core
{

// generic fallthrough
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<auto type, typename Src, typename Dst, typename GridLayout>
void ParticlesAppender<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    Src const& src, Dst& dst, GridLayout const& /*layout*/)
{
    using enum LayoutMode;
    static_assert(src_layout_mde == dst_layout_mde);
    static_assert(src_alloc_mde == dst_alloc_mde);
    static_assert(any_in(src_layout_mde, AoS, SoA));

    auto const oldSize = dst.size();
    dst.resize(oldSize + src.size());

    if constexpr (src_layout_mde == AoS)
    {
        mem::copy<src_alloc_mde>(dst.data() + oldSize, src.data(), src.size());
    }
    else if constexpr (src_layout_mde == SoA)
    {
        auto dst_tuple = dst.as_tuple();
        auto src_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(dst_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<src_alloc_mde>(a.data() + oldSize, b, src.size());
        });
    }
    else
        throw std::runtime_error("no");

    // dst.resize(dst.size() + src.size());

    // ParticleArrayService::sync<2, type>(dst);
}




} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING */
