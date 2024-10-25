#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include <type_traits>

namespace PHARE::core::detail
{

// template<std::size_t dim, auto lm, auto sm, auto am, std::uint8_t impl>
// struct ParticleArraySelectorImpl;


template<typename ParticleArray_t, std::uint8_t _impl>
class ParticleArraySelector
{
    auto static constexpr dimension    = ParticleArray_t::dim;
    auto static constexpr alloc_mode   = ParticleArray_t::alloc_mode;
    auto static constexpr layout_mode  = ParticleArray_t::layout_mode;
    auto static constexpr storage_mode = ParticleArray_t::storage_mode;

public:
    // using dispatch
    //     = ParticleArraySelectorImpl<dimension, layout_mode, storage_mode, alloc_mode, _impl>;
    using box_t = Box<int, dimension>;

    template<typename ParticleArray_t0>
    void operator()(ParticleArray_t0& dst, box_t const& select) const; /*
   {
       if constexpr (ParticleArray_t::is_mapped)
           from.export_particles(select, to);
   }*/

    template<typename ParticleArray_t0, typename Transformer>
    void operator()(ParticleArray_t0& dst, box_t const& select, Transformer&& fn) const;

    // box_t const domain_box;
    ParticleArray_t const& from;
    // ParticleArray& to;
    std::size_t start = 0, end = from.size();
};

template<auto am, auto lm, std::uint8_t impl, typename ParticleArray_t>
struct SelectorImpl;


template<typename ParticleArray_t>
struct SelectorImpl<AllocatorMode::CPU, LayoutMode::AoS, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
        // PHARE_LOG_LINE_SS("");
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct SelectorImpl<AllocatorMode::CPU, LayoutMode::SoA, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
        // PHARE_LOG_LINE_SS("");
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct SelectorImpl<AllocatorMode::GPU_UNIFIED, LayoutMode::AoS, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
        // PHARE_LOG_LINE_SS("");
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};

template<typename ParticleArray_t>
struct SelectorImpl<AllocatorMode::GPU_UNIFIED, LayoutMode::SoA, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&)
    {
        // PHARE_LOG_LINE_SS("");
    }

    template<typename ParticleArray_t0, typename Transformer>
    static void select(ParticleArray_t const&, ParticleArray_t0&, box_t const&, Transformer&&)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
    }
};


template<typename ParticleArray_t>
struct SelectorImpl<AllocatorMode::CPU, LayoutMode::AoSPC, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
        // PHARE_LOG_LINE_SS("");
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
struct SelectorImpl<AllocatorMode::CPU, LayoutMode::SoAPC, /*impl = */ 0, ParticleArray_t>
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
struct SelectorImpl<AllocatorMode::GPU_UNIFIED, LayoutMode::AoSPC, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);
        // PHARE_LOG_LINE_SS("");
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
struct SelectorImpl<AllocatorMode::GPU_UNIFIED, LayoutMode::SoAPC, /*impl = */ 0, ParticleArray_t>
{
public:
    using box_t = Box<int, ParticleArray_t::dimension>;

    template<typename ParticleArray_t0>
    static void select(ParticleArray_t const& src, ParticleArray_t0& dst, box_t const& select)
    {
        static_assert(std::is_same_v<ParticleArray_t, ParticleArray_t0>);

        // PHARE_LOG_LINE_SS("");

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

template<typename ParticleArray_t, std::uint8_t impl>
template<typename ParticleArray_t0>
void ParticleArraySelector<ParticleArray_t, impl>::operator()(ParticleArray_t0& dst,
                                                              box_t const& select) const
{
    SelectorImpl<ParticleArray_t::alloc_mode, ParticleArray_t::layout_mode, impl,
                 ParticleArray_t>::select(from, dst, select);
}


template<typename ParticleArray_t, std::uint8_t impl>
template<typename ParticleArray_t0, typename Transformer>
void ParticleArraySelector<ParticleArray_t, impl>::operator()(ParticleArray_t0& dst,
                                                              box_t const& select,
                                                              Transformer&& fn) const
{
    SelectorImpl<ParticleArray_t::alloc_mode, ParticleArray_t::layout_mode, impl,
                 ParticleArray_t>::select(from, dst, select, std::forward<Transformer>(fn));
}


} // namespace PHARE::core::detail

// #include "particles_selecting_cpu.hpp"

// #if PHARE_HAVE_MKN_GPU
// #include "particles_selecting_gpu_mkn.hpp"
// #else // no impl (yet)
// namespace PHARE::core::detail
// {
// template<typename ParticleArray_t>
// class ParticleArraySelector<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray_t>;
// // template<typename ParticleArray_t>
// // class ParticleSelector<AllocatorMode::GPU, /*impl = */ 0, ParticleArray_t>;
// } // namespace PHARE::core::detail
// #endif

#endif /*PHARE_CORE_DATA_PARTICLES_SELECTING_PARTICLES_SELECTING_HPP*/
