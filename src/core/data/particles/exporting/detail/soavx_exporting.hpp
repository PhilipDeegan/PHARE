#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOAVX_EXPORTING
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOAVX_EXPORTING

// SoAVX: untiled and untracked, particles leaving a box are partitioned over the zipped columns
//  (thrust), then moved to dst / dropped from the tail. Needs thrust, as SoAVX iteration does
//  GPU_UNIFIED partitions and moves on the device (thrust::device, device copies)

#include "core/def.hpp"
#include "core/def/phare_config.hpp"
#include "core/vector.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/memory.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>
#include <iterator>
#include <stdexcept>

#if PHARE_HAVE_THRUST
#include "core/data/particles/arrays/particle_array_soa_thrust.hpp"
#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#endif

namespace PHARE::core
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesExporter;

template<auto alloc_mde>
struct ParticlesExporter<LayoutMode::SoAVX, alloc_mde>
{
    auto constexpr static layout_mode = LayoutMode::SoAVX;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Src, std::size_t dim>
    void delete_particles_in(Src& src, Box<int, dim> const& box)
    {
        src.resize(partition_(src, [=] _PHARE_ALL_FN_(auto const& cell) { return !isIn(cell, box); }));
    }

    template<typename Src, typename Boxes>
    void delete_particles_in(Src& src, Boxes const& boxes)
    {
        auto const list  = box_list_(boxes);
        auto const view  = list.view();
        src.resize(partition_(src, [=] _PHARE_ALL_FN_(auto const& cell) { return !view.contains(cell); }));
    }

    template<typename Src, std::size_t dim>
    void delete_particles_not_in(Src& src, Box<int, dim> const& box)
    {
        src.resize(partition_(src, [=] _PHARE_ALL_FN_(auto const& cell) { return isIn(cell, box); }));
    }

    template<typename Src, typename Boxes>
    void delete_particles_not_in(Src& src, Boxes const& boxes)
    {
        auto const list = box_list_(boxes);
        auto const view = list.view();
        src.resize(partition_(src, [=] _PHARE_ALL_FN_(auto const& cell) { return view.contains(cell); }));
    }

    template<typename Dst, typename Src, std::size_t dim>
    void move_in_domain(Dst& dst, Src& src, Box<int, dim> const& domain_box)
    {
        auto const kept = partition_(
            src, [=] _PHARE_ALL_FN_(auto const& cell) { return !isIn(cell, domain_box); });
        move_tail_(dst, src, kept);
    }

    template<typename Dst, typename Src, std::size_t dim>
    void move_in_ghost_layer(Dst& dst, Src& src, Box<int, dim> const& domain_box,
                             Box<int, dim> const& ghost_box)
    {
        auto const list = box_list_(ghost_box.remove(domain_box));
        auto const view = list.view();
        auto const kept = partition_(
            src, [=] _PHARE_ALL_FN_(auto const& cell) { return !view.contains(cell); });
        move_tail_(dst, src, kept);
    }

    template<typename Dst, typename Src, std::size_t dim, typename Boxes>
    void move_in_ghost_layer(Dst& dst, Src& src, Box<int, dim> const& domain_box,
                             Boxes const& ghost_boxes)
    {
        for (auto const& gb : ghost_boxes)
            move_in_ghost_layer(dst, src, domain_box, gb);
    }

private:
    // a box list the partition predicate can read where it runs: copied into this allocator's
    //  memory (unified for GPU_UNIFIED), read through pointer and size
    template<typename Box_t>
    struct BoxListView
    {
        Box_t const* boxes;
        std::size_t size;

        bool contains(auto const& cell) const _PHARE_ALL_FN_
        {
            for (std::size_t i = 0; i < size; ++i)
                if (isIn(cell, boxes[i]))
                    return true;
            return false;
        }
    };

    template<typename Boxes>
    static auto box_list_(Boxes const& boxes)
    {
        using Box_t    = std::decay_t<typename Boxes::value_type>;
        using Vector_t = typename PHARE::Vector<Box_t, alloc_mode>::vector_t;
        struct List
        {
            Vector_t boxes;
            auto view() const { return BoxListView<Box_t>{boxes.data(), boxes.size()}; }
        };
        List list{PHARE::Vector<Box_t, alloc_mode>::make(0)};
        list.boxes.reserve(boxes.size());
        for (auto const& box : boxes)
            list.boxes.push_back(box);
        return list;
    }

    // particles for which keep(iCell) first: returns how many, on the device for GPU_UNIFIED
    static std::size_t partition_([[maybe_unused]] auto& src, [[maybe_unused]] auto const keep)
    {
#if PHARE_HAVE_THRUST
        using Iterator_t = SoAVXIteratorAdaptor;
        auto const it    = Iterator_t::make(src);
        auto const pred  = [=] _PHARE_ALL_FN_(auto&& z) { return keep(Iterator_t::iCell(z)); };
        auto const pivot = [&]() {
            if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
                return thrust::partition(thrust::device, it, it + src.size(), pred);
            else
                return thrust::partition(thrust::host, it, it + src.size(), pred);
        }();
        return static_cast<std::size_t>(std::distance(it, pivot));
#else
        throw std::runtime_error("SoAVX partitioning needs thrust");
#endif
    }

    // src's particles from `from` onwards appended to dst, src keeps [0, from): dst presized,
    //  then per column the tail copied at its offset (on the device for GPU_UNIFIED)
    static void move_tail_(auto& dst, auto& src, std::size_t const from)
    {
        std::size_t const n      = src.size() - from;
        std::size_t const offset = dst.size();
        if (n > 0)
        {
            auto dst_columns       = dst.as_tuple();
            auto const src_columns = src.as_tuple();
            for_N<std::tuple_size_v<decltype(dst_columns)>>([&](auto c) {
                auto& d       = std::get<c>(dst_columns);
                auto const& s = std::get<c>(src_columns);
                PHARE::resize(d, offset + n);
                mem::copy<alloc_mode>(d.data() + offset, s.data() + from, n);
            });
        }
        src.resize(from);
    }
};

// SoAVXChunked: the same column moves, every one ends in a resize, which unsorts: the offsets
//  are rebuilt by the following on_moved / on_appended
template<auto alloc_mde>
struct ParticlesExporter<LayoutMode::SoAVXChunked, alloc_mde>
    : public ParticlesExporter<LayoutMode::SoAVX, alloc_mde>
{
    auto constexpr static layout_mode = LayoutMode::SoAVXChunked;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOAVX_EXPORTING */
