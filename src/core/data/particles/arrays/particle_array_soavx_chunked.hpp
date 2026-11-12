#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_CHUNKED_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_CHUNKED_HPP

// SoAVXChunked: SoAVX sorted by cell, with per cell offsets (CSR), as AoSMapped is AoS + CellMap
//
//  offsets_[c] is the first particle of cell c, offsets_[c + 1] - offsets_[c] its count, cells c
//  flattened in C order (last direction fastest) over shell() = grow(box, 1) (as AoSMapped's
//  CellMap, the safe outer shell). A span of cells along the last direction is then one
//  contiguous range of particles: selecting / exporting a box is a list of ranges.
//
//  The offsets are rebuilt by sort() only: a counting sort (histogram with a slot per particle,
//  exclusive scan, one scatter of all columns into a second set of columns, swapped in).
//  Any mutation (push, append, erase, resize) leaves them stale: sorted() is false until the
//  next sort(), range accessors assert sorted(). GPU_UNIFIED sorts on the device (thrust), the
//  order of particles within a cell is then not deterministic (atomic slots), CPU is stable.

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/operators.hpp"
#include "core/utilities/memory.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"

#include <tuple>
#include <cstdint>
#include <numeric>

#if PHARE_HAVE_THRUST
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#endif

namespace PHARE::core
{

// per cell offsets type: particles per patch fit in 32 bits, and device atomics on it are cheap
using chunk_offset_t = std::uint32_t;


// the offsets view shared by the vector and the span: shell box + offsets pointer
template<std::size_t dim>
struct SoAVXChunks
{
    using box_t = Box<int, dim>;

    // flat C order index of cell (x[, y[, z]]) in shell
    static std::size_t flat(box_t const& shell, auto const& cell) _PHARE_ALL_FN_
    {
        std::size_t idx = 0;
        for (std::size_t d = 0; d < dim; ++d)
            idx = idx * (shell.upper[d] - shell.lower[d] + 1) + (cell[d] - shell.lower[d]);
        return idx;
    }

    // calls fn(start, size) per contiguous particle range of the cells of box (clipped to
    //  shell), adjacent ranges merged, empty ones skipped
    static void for_each_range(box_t const& shell, chunk_offset_t const* offsets,
                               box_t const& box, auto&& fn)
    {
        auto const clipped = box * shell;
        if (!clipped)
            return;

        std::size_t start = 0, size = 0;
        for_each_span(*clipped, [&](auto const& span_start, std::uint32_t const span_size) {
            auto const c0 = flat(shell, span_start);
            auto const s  = offsets[c0];
            auto const n  = offsets[c0 + span_size] - s;
            if (n == 0)
                return;
            if (size > 0 and start + size == s)
                size += n;
            else
            {
                if (size > 0)
                    fn(start, size);
                start = s;
                size  = n;
            }
        });
        if (size > 0)
            fn(start, size);
    }
};


template<std::size_t dim, auto alloc_mode_, auto _const_ = 0>
struct SoAVXChunkedSpan : public SoAVXSpan<dim, alloc_mode_, _const_>
{
    using Super = SoAVXSpan<dim, alloc_mode_, _const_>;
    using box_t = Box<int, dim>;

    SoAVXChunkedSpan() = default;

    template<typename Particles_t>
    SoAVXChunkedSpan(Particles_t& particles)
        : Super{particles, 0, particles.size()}
        , box_{particles.box()}
        , offsets_{particles.offsets_data()}
        , sorted_{particles.sorted()}
    {
    }

    auto& box() const _PHARE_ALL_FN_ { return box_; }
    auto offsets_data() const _PHARE_ALL_FN_ { return offsets_; }
    bool sorted() const _PHARE_ALL_FN_ { return sorted_; }

    box_t box_;
    chunk_offset_t const* offsets_ = nullptr;
    bool sorted_                   = false;
};


template<std::size_t dim, auto alloc_mode_>
struct SoAVXChunkedVector : public SoAVXVector<dim, alloc_mode_>
{
    using Super = SoAVXVector<dim, alloc_mode_>;
    using box_t = Box<int, dim>;

    template<typename T>
    using container_t = Super::template container_t<T>;

    SoAVXChunkedVector(box_t const& box = {}, std::size_t const size = 0)
        : Super(size)
        , box_{box}
    {
    }

    auto& box() const { return box_; }
    auto offsets_data() const { return offsets_.data(); }
    bool sorted() const { return sorted_ and sorted_size_ == this->size(); }

    // the offsets are stale after any change to the particles, until the next sort()
    void unsort() { sorted_ = false; }

    box_t box_;
    container_t<chunk_offset_t> offsets_; // ncells + 1, see top of file
    container_t<chunk_offset_t> slots_;   // per particle slot in its cell, sort() scratch
    Super scratch_;                       // the columns sort() scatters into, then swaps in
    bool sorted_              = false;
    std::size_t sorted_size_ = 0;
};


template<typename Super_>
class SoAVXChunkedParticles : public SoAVXParticles<Super_>
{
    using Super  = SoAVXParticles<Super_>;
    using Chunks = SoAVXChunks<Super::dimension>;

public:
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr layout_mode  = LayoutMode::SoAVXChunked;
    auto static constexpr storage_mode = Super::storage_mode;

    using box_t = Box<int, dimension>;
    using Super::size;

    template<typename... Args>
    SoAVXChunkedParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    SoAVXChunkedParticles(SoAVXChunkedParticles const&)            = default;
    SoAVXChunkedParticles(SoAVXChunkedParticles&&)                 = default;
    SoAVXChunkedParticles& operator=(SoAVXChunkedParticles const&) = default;
    SoAVXChunkedParticles& operator=(SoAVXChunkedParticles&&)      = default;

    // the cells the offsets cover
    auto shell() const { return grow(this->box(), 1); }

    // particles of cell, as {start, size}
    auto cell_range(auto const& cell) const
    {
        PHARE_ASSERT(this->sorted());
        auto const* offsets = this->offsets_data();
        auto const c        = Chunks::flat(shell(), cell);
        return std::pair<std::size_t, std::size_t>{offsets[c], offsets[c + 1] - offsets[c]};
    }

    // fn(start, size) per contiguous range of the particles in box, see SoAVXChunks
    void for_each_range(box_t const& box, auto&& fn) const
    {
        PHARE_ASSERT(this->sorted());
        Chunks::for_each_range(shell(), this->offsets_data(), box, fn);
    }

    std::size_t count(box_t const& box) const
    {
        std::size_t n = 0;
        for_each_range(box, [&](auto, auto const size) { n += size; });
        return n;
    }

    // counting sort by cell over shell(), all columns moved together, offsets rebuilt
    void sort()
        requires(storage_mode == StorageMode::VECTOR);

    void resize(std::size_t const s)
        requires(storage_mode == StorageMode::VECTOR)
    {
        this->unsort();
        Super_::resize(s);
    }

    template<typename Particle_t>
    void push_back(Particle_t const& particle)
        requires(storage_mode == StorageMode::VECTOR)
    {
        this->unsort();
        Super::push_back(particle);
    }

    void clear()
        requires(storage_mode == StorageMode::VECTOR)
    {
        this->unsort();
        Super_::clear();
    }

    // one particle: anything with weight()/iCell()/delta()/v(), or Particle_t's arguments
    template<typename... Args>
    void emplace_back(Args&&... args)
        requires(storage_mode == StorageMode::VECTOR)
    {
        if constexpr (sizeof...(Args) == 1)
            push_back(args...);
        else
            push_back(Particle_t{args...});
    }

    void reserve(std::size_t const s)
        requires(storage_mode == StorageMode::VECTOR)
    {
        std::apply([&](auto&... column) { ((column.reserve(s)), ...); }, this->as_tuple());
    }

    // the offsets are kept eagerly, as AoSMapped's CellMap: rebuilt after every change
    //  through these hooks (appends, the post push sync), queries assert sorted()
    template<auto type>
    void on_appended()
        requires(storage_mode == StorageMode::VECTOR)
    {
        sort();
    }

    template<auto type>
    void on_moved(auto&&... /*stream*/)
        requires(storage_mode == StorageMode::VECTOR)
    {
        sort();
    }

    // that's particles appended, column by column: one bulk copy per column whatever that's
    //  allocator (e.g. a CPU staging array into GPU_UNIFIED). Unsorted until on_appended
    template<typename That>
    void append_columns(That const& that)
        requires(storage_mode == StorageMode::VECTOR)
    {
        this->unsort();
        std::size_t const offset = size();
        std::size_t const n      = that.size();
        if (n == 0)
            return;

        auto columns            = this->as_tuple();
        auto const that_columns = that.as_tuple();
        for_N<std::tuple_size_v<decltype(columns)>>([&](auto c) {
            auto& column = std::get<c>(columns);
            PHARE::resize(column, offset + n);
            mem::copy<alloc_mode>(column.data() + offset, std::get<c>(that_columns).data(), n);
        });
    }

    std::size_t nbr_particles_in(box_t const& box) const { return count(box); }

    std::size_t nbr_particles_in(std::array<int, dimension> const& cell) const
    {
        return cell_range(cell).second;
    }

private:
    using Particle_t = Particle<dimension>;
};


template<typename Super_>
void SoAVXChunkedParticles<Super_>::sort()
    requires(storage_mode == StorageMode::VECTOR)
{
    PHARE_LOG_SCOPE(2, "SoAVXChunked::sort");
    auto static constexpr gpu = alloc_mode == AllocatorMode::GPU_UNIFIED;
    using Op                  = Operators<chunk_offset_t, /*atomic=*/gpu, /*GPU=*/gpu>;

    auto const shell        = this->shell();
    std::size_t const n     = size();
    std::size_t const cells = shell.size();

    PHARE::resize(this->offsets_, cells + 1, /*copy=*/false);
    PHARE::resize(this->slots_, n, /*copy=*/false);
    PHARE::resize(this->scratch_.weight_, n, false);
    for_N<dimension>([&](auto d) {
        PHARE::resize(this->scratch_.iCell_[d], n, false);
        PHARE::resize(this->scratch_.delta_[d], n, false);
    });
    for_N<3>([&](auto d) { PHARE::resize(this->scratch_.v_[d], n, false); });

    auto* const offsets = this->offsets_.data();
    auto* const slots   = this->slots_.data();
    auto cells_of       = for_N<dimension, for_N_R_mode::make_array>(
        [&](auto d) -> int const* { return this->iCell_[d].data(); });
    auto const cell_idx = [=] _PHARE_ALL_FN_(std::size_t const i) {
        auto const cell = for_N<dimension, for_N_R_mode::make_array>(
            [&](auto d) { return cells_of[d][i]; });
        PHARE_ASSERT(isIn(Point{cell}, shell));
        return Chunks::flat(shell, cell);
    };

    // particles' columns in, scratch columns out, flat lists of pointers
    auto const columns = [](auto& v) {
        return std::apply([](auto&... c) { return std::make_tuple(c.data()...); }, v.as_tuple());
    };
    auto const src = columns(*this);
    auto const dst = columns(this->scratch_);

    auto const count = [=] _PHARE_ALL_FN_(std::size_t const i) {
        slots[i] = Op{offsets[cell_idx(i)]}.increment_return_old();
    };
    auto const scatter = [=] _PHARE_ALL_FN_(std::size_t const i) {
        auto const to = offsets[cell_idx(i)] + slots[i];
        for_N<std::tuple_size_v<decltype(src)>>(
            [&](auto c) { std::get<c>(dst)[to] = std::get<c>(src)[i]; });
    };

    if constexpr (gpu)
    {
#if PHARE_HAVE_THRUST
        auto const first = thrust::counting_iterator<std::size_t>(0);
        thrust::fill(thrust::device, offsets, offsets + cells + 1, chunk_offset_t{0});
        thrust::for_each(thrust::device, first, first + n, count);
        thrust::exclusive_scan(thrust::device, offsets, offsets + cells + 1, offsets);
        thrust::for_each(thrust::device, first, first + n, scatter);
#else
        throw std::runtime_error("SoAVXChunked GPU sort needs thrust");
#endif
    }
    else
    {
        std::fill(offsets, offsets + cells + 1, chunk_offset_t{0});
        for (std::size_t i = 0; i < n; ++i)
            count(i);
        std::exclusive_scan(offsets, offsets + cells + 1, offsets, chunk_offset_t{0});
        for (std::size_t i = 0; i < n; ++i)
            scatter(i);
    }

    auto mine   = this->as_tuple();
    auto theirs = this->scratch_.as_tuple();
    for_N<std::tuple_size_v<decltype(mine)>>(
        [&](auto c) { std::get<c>(mine).swap(std::get<c>(theirs)); });

    this->sorted_      = true;
    this->sorted_size_ = n;
}


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using SoAVXChunkedVectorParticles = SoAVXChunkedParticles<SoAVXChunkedVector<dim, alloc_mode>>;


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_CHUNKED_HPP */
