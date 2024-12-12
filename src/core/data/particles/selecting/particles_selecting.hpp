#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP

#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <type_traits>

namespace PHARE::core::detail
{



template<auto am, auto lm, std::uint8_t impl, typename ParticleArray_t>
struct ParticleArraySelector;


template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::CPU, LayoutMode::AoS, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::CPU, LayoutMode::SoA, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::GPU_UNIFIED, LayoutMode::AoS, /*impl = */ 0,
                             ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::GPU_UNIFIED, LayoutMode::SoA, /*impl = */ 0,
                             ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};


template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::CPU, LayoutMode::AoSPC, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);

        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select);
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

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select,
                       Transformer&& fn)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select - fn); // BEWARE
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
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::CPU, LayoutMode::SoAPC, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);

        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select);
        assert(lcl_src_box.shape() == lcl_dst_box.shape());
        auto src_it = lcl_src_box.begin();
        auto dst_it = lcl_dst_box.begin();
        for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
        {
            auto& sv = src(*src_it);
            auto& dv = dst(*dst_it);
            dv.reserve(dv.size() + sv.size());
            dv.emplace_back(sv);
        }
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select,
                       Transformer&& fn)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);

        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select - fn); // BEWARE
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
                p.iCell() = (Point{p.iCell()} + fn).toArray();
                dv.emplace_back(p);
            }
        }
    }
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::GPU_UNIFIED, LayoutMode::AoSPC, /*impl = */ 0,
                             ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);

        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select);
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

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct ParticleArraySelector<AllocatorMode::GPU_UNIFIED, LayoutMode::SoAPC, /*impl = */ 0,
                             ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);



        auto const lcl_src_box = src.local_box(select);
        auto const lcl_dst_box = dst.local_box(select);
        assert(lcl_src_box.shape() == lcl_dst_box.shape());
        auto src_it = lcl_src_box.begin();
        auto dst_it = lcl_dst_box.begin();
        for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
        {
            auto& sv = src(*src_it);
            auto& dv = dst(*dst_it);
            dv.reserve(dv.size() + sv.size());
            dv.emplace_back(sv);
        }
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};


} // namespace PHARE::core::detail


#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP*/
