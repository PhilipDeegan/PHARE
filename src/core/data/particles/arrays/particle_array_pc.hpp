#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP

#include "core/data/particles/particle.hpp"
#include "core/operators.hpp"
#include "core/utilities/span.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "particle_array_aos.hpp"

#include <iterator>
#include <unordered_set>

namespace PHARE::core
{

template<typename Particles, std::uint8_t impl_ = 0>
class PerCellSpan
{
protected:
    using SIZE_T = unsigned long long int;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;
    auto static constexpr impl_v     = impl_;
    using This                       = PerCellSpan<Particles, impl_v>;
    using lobox_t                    = Box<std::uint32_t, dim>;
    using per_cell_particles         = Particles;

private:
    using locell_t = std::array<std::uint32_t, dim>;

    template<typename PerCellArray>
    auto resolve(PerCellArray& arr)
    {
        if constexpr (PerCellArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
        {
            arr.check();
            return *arr.particles_views_;
        }
    }

    template<typename PerCellArray>
    auto resolve_gaps(PerCellArray& arr)
    {
        if constexpr (PerCellArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
        {
            arr.check();
            return *arr.gap_views_;
        }
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = Particle<dim>;

    template<typename PerCellArray>
    PerCellSpan(PerCellArray& arr)
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , off_sets_{arr.off_sets_}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cap_{arr.cap_}
        , left_{arr.left_}
        , p2c_{arr.p2c_.data(), arr.p2c_.size()}
        , size_{arr.size()}
        , box_{arr.box_}
        , ghost_box_{arr.ghost_box()}
        , safe_box_{arr.safe_box_}
        , local_ghost_box_{arr.local_box()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    auto size(std::size_t const& idx) const _PHARE_ALL_FN_ { return particles_.data()[idx].size(); }

    // void resize(std::size_t s) { particles_.s = s; }

    auto& box() const _PHARE_ALL_FN_ { return box_; }
    auto& ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }
    auto local_cell(std::array<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return as_local_cell(safe_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return local_cell(icell.toArray());
    }

    auto& local_box() const _PHARE_ALL_FN_ { return local_ghost_box_; }

    auto local_box(Box<int, dim> const& from) const _PHARE_ALL_FN_
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) _PHARE_ALL_FN_ { return particles_(cell); }
    auto& operator()(locell_t const& cell) const _PHARE_ALL_FN_ { return particles_(cell); }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain, typename... Args>
    void sync(Args&&... args) _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_add_new() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left() _PHARE_ALL_FN_;

    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        size_ = 0;
    }

    void reset_p2c(std::array<std::uint32_t, dim>* cells, std::size_t size)
    {
        p2c_ = Span<std::array<std::uint32_t, dim>>{cells, size};
    }


protected:
    NdArrayView<dim, Particles> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> off_sets_, gap_idx_, add_into_, cap_, left_;
    Span<std::array<std::uint32_t, dim>> p2c_;
    std::size_t size_;

    Box<int, dim> box_, ghost_box_, safe_box_;
    lobox_t local_ghost_box_;


}; // PerCellSpan


template<typename Particles, std::uint8_t impl_ = 0>
class PerCellVector
{
    using This   = PerCellVector<Particles, impl_>;
    using SIZE_T = unsigned long long int; // cuda issues

    template<typename P, std::uint8_t i>
    friend class PerCellSpan;

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr impl_v       = impl_;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr layout_mode  = Particles::layout_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = Particle<dim>;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::Span_t;
    using per_cell_particles           = Particles;

    template<typename T>
    using vec_helper = PHARE::Vector<T, alloc_mode, 1>;

    using size_t_vec_helper   = vec_helper<std::size_t>;
    using size_t_vector       = typename size_t_vec_helper::vector_t;
    using particle_vec_helper = vec_helper<Particle_t>;
    using particle_vector     = typename particle_vec_helper::vector_t;

    PerCellVector(box_t const& box = {}, std::size_t ghost_cells = 0)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
        , safe_box_{grow(box, ghost_cells + 1)}
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        left_.zero();
        cap_.zero();
    }

    PerCellVector(PerCellVector const& from)            = default;
    PerCellVector(PerCellVector&& from)                 = default;
    PerCellVector& operator=(PerCellVector&& from)      = default;
    PerCellVector& operator=(PerCellVector const& from) = default;

    template<auto type = ParticleType::Domain>
    auto& reserve_ppc(std::size_t const& ppc);

    auto size() const { return total_size; }
    auto size(std::array<std::uint32_t, dim> const& icell) const { return cell_size_(icell); }
    auto size(std::size_t const& idx) const { return cell_size_.data()[idx]; }

    template<typename Iterator>
    auto erase(Iterator first, Iterator last)
    {
        return particles_.erase(particles_.begin() + first.curr_pos,
                                particles_.begin() + last.curr_pos);
    }

    void _inc(locell_t const& locell)
    {
        ++total_size;
        ++cell_size_(locell);
    }

    template<bool inc_ = true>
    void /*Particle_t&*/ emplace_back(Particle_t const& p)
    {
        // auto& np = get_vec(particles_(local_cell(p.iCell()))).emplace_back(p);
        auto const locell = local_cell(p.iCell());
        /*auto& np = */ particles_(locell).emplace_back(p);
        if constexpr (inc_)
            _inc(locell);
        // return np;
    }

    template<bool inc_ = true>
    void emplace_back(Particles& dst, Particles const& src, std::size_t const& idx)
    {
        dst.emplace_back(src, idx);
    }

    template<typename V>
    static auto& get_vec(V& v)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(return mkn::gpu::as_super(v));
        }
        else
        {
            return v;
        }
    }

    void push_back(Particle_t&& p) { emplace_back(p); }
    void push_back(Particle_t const& p) { emplace_back(p); }

    void reset_views()
    {
        particles_views_ = generate_from<alloc_mode>(
            [&](auto const i) {
                return PSpan_t{*(particles_.data() + i), *(cell_size_.data() + i)};
            },
            particles_);
        gap_views_ = generate_from<alloc_mode>(
            [&](auto const i) { return make_span(*(gaps_.data() + i)); }, gaps_);
    }

    void reset_index_wrapper_map();

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }
    auto& safe_box() const { return safe_box_; } // allocating box


    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(safe_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }
    auto local_box() const
    {
        return box_from_zero_to_upper_minus_one(
            safe_box_.shape().template toArray<std::uint32_t>());
    }
    auto local_box(Box<int, dim> const& from) const
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync(); // after move

    template<auto type>
    void trim();



    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        total_size = 0;
    }


    auto& insert(PerCellVector const& src);

    auto& insert_domain_from(PerCellVector const& src);


    void replace_from(This const& that)
    {
        throw std::runtime_error("fix");
        if (this == &that) // just in case
            return;

        cell_size_ = that.cell_size_;
        particles_ = that.particles_;
        box_       = that.box_;
        ghost_box_ = that.ghost_box_;
        p2c_       = that.p2c_;
    }

    void cap()
    {
        PHARE_LOG_SCOPE(1, "PerCellVector::cap");
        if constexpr (any_in(impl_v, 1, 2))
            for (auto& p : particles_)
                resize(p, p.capacity());
    }



    void static resize(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { resize(v, s, copy); });
        else
            resize(ps.particles_, s, copy);
    }

    template<typename V>
    void static resize(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::resize(v, s, copy));

        else
            v.resize(s);
    }

    void static reserve(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { reserve(v, s, copy); });
        else
            reserve(ps.particles_, s, copy);
    }

    template<typename V>
    void static reserve(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::reserve(v, s, copy));

        else
            v.reserve(s);
    }

    template<typename View_t>
    void reset_p2c(View_t& view)
    {
        view.reset_p2c(p2c_.data(), p2c_.size());
    }

protected:
    bool static constexpr c_order = true;

    std::size_t ghost_cells_;

    Box<int, dim> box_, ghost_box_, safe_box_;

    NdArrayVector<dim, Particles, c_order, alloc_mode> particles_{local_box().shape()};

    // only used for GPU

    NdArrayVector<dim, PSpan_t, c_order, alloc_mode> particles_views_{local_box().shape()};

    typename vec_helper<std::array<std::uint32_t, dim>>::vector_t p2c_;
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> off_sets_{local_box().shape()};

    NdArrayVector<dim, size_t_vector, c_order, alloc_mode> gaps_{local_box().shape()};
    NdArrayVector<dim, Span<std::size_t>, c_order, alloc_mode> gap_views_{local_box().shape()};

    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> gap_idx_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> add_into_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> left_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> cap_{local_box().shape()};
    NdArrayVector<dim, std::size_t, c_order, alloc_mode> cell_size_{local_box().shape()};

    std::size_t total_size = 0;



    template<auto type>
    void sync_cpu_gaps_and_tmp();
    template<auto type>
    void sync_gpu_gaps_and_tmp();

    template<auto type>
    void sync_gpu_gaps_and_tmp_impl0();

    template<auto type>
    void sync_gpu_gaps_and_tmp_impl1();

    static void on_box(auto&& box, auto&& fn)
    {
        for (auto const& bix : box)
            fn(bix);
    };

    static void on_box_list(auto&& boxlist, auto&& fn)
    {
        for (auto const& box : boxlist)
            on_box(box, fn);
    };

    void on_domain(auto&& fn) const { on_box(local_box(box()), fn); };
    void on_ghost_box(auto&& fn) { on_box(local_box(), fn); };
    void on_ghost_layer(auto&& fn) const { on_box_list(local_box().remove(local_box(box())), fn); };
    void on_ghost_layer_plus_1_domain(auto&& fn) const
    {
        on_box_list(local_box().remove(shrink(local_box(box()), 1)), fn);
    };
    void on_ghost_layer_plus_2_domain(auto&& fn) const
    {
        on_box_list(local_box().remove(shrink(local_box(box()), 2)), fn);
    };

}; // PerCellVector<Particles>



template<typename Particles, std::uint8_t impl>
auto& PerCellVector<Particles, impl>::insert(PerCellVector const& src)
{
    std::size_t added = 0;
    for (auto const& bix : local_box(box()))
    {
        auto& from = src(bix);
        auto& to   = (*this)(bix);
        added += from.size();
        to.reserve(to.size() + from.size());
        std::copy(from.begin(), from.end(), std::back_inserter(to));
    }
    if (added)
        sync<2>();

    return *this;
}

template<typename Particles, std::uint8_t impl>
auto& PerCellVector<Particles, impl>::insert_domain_from(PerCellVector const& src)
{
    auto const on_box = [&](auto&& box, auto&& fn) {
        for (auto const& bix : box)
            fn(bix);
    };
    auto const on_box_list = [&](auto&& boxlist, auto&& fn) {
        for (auto const& box : boxlist)
            on_box(box, fn);
    };

    std::size_t added = 0;

    // this->check();
    // src.check();

    on_box_list(local_box(box()).remove(shrink(local_box(box()), 1)), [&](auto const& bix) {
        auto& from = src(bix);
        auto& to   = (*this)(bix);
        added += from.size();
        to.reserve(to.size() + from.size());

        to.emplace_back(from);
        // for (auto const& p : from)
        //     to.emplace_back(p);

        // std::copy(from.begin(), from.end(), std::back_inserter(to));
    });

    // this->check();


    if (added)
        sync<2>();

    // this->check();

    return *this;
}

template<typename Particles, std::uint8_t impl>
template<auto type>
auto& PerCellVector<Particles, impl>::reserve_ppc(std::size_t const& ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    std::size_t const additional = ppc < 50 ? 10 : ppc * .2; // 20% overallocate
    std::size_t const buffered   = ppc + additional;

    if constexpr (type == ParticleType::Domain)
    {
        on_domain([&](auto const& bix) { particles_(bix).reserve(buffered); });
        on_ghost_box([&](auto const& bix) { reserve(gaps_(bix), additional); });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(additional); });
    }

    if constexpr (type == ParticleType::Ghost)
    {
        on_ghost_layer_plus_2_domain([&](auto const& bix) {
            particles_(bix).reserve(additional);
            reserve(gaps_(bix), additional);
        });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(buffered); });
    }

    return *this;
}



template<typename Particles, std::uint8_t impl>
template<auto type>
void PerCellVector<Particles, impl>::trim()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    if constexpr (type == ParticleType::Ghost)
        for (auto const& ghost_layer_box : local_box().remove(local_box(box_)))
            for (auto const& bix : ghost_layer_box)
            {
                cell_size_(bix) = 0;
                particles_(bix).clear();
                gaps_(bix).clear();
            }

    if constexpr (type == ParticleType::Domain)
    {
        throw std::runtime_error("NO");
    }
}

template<typename Particles, std::uint8_t impl>
template<auto type>
void PerCellVector<Particles, impl>::sync_cpu_gaps_and_tmp()
{
    // PHARE_LOG_LINE_STR("sync_cpu_gaps_and_tmp " << magic_enum::enum_name(type));
    auto const lbox = local_box();

    for (auto const& bix : lbox)
    {
        auto& real = particles_(bix);
        auto& gaps = gaps_(bix);
        int diff   = real.size() - cell_size_(bix); // new additions
        if (diff < 0)
        {
            PHARE_LOG_LINE_STR(real.size()
                               << " " << gaps.size() << " " << cell_size_(bix) << " " << diff);
        }
        PHARE_ASSERT(diff >= 0);
        while (gaps.size() and diff)
        {
            PHARE_ASSERT(gaps.back() < real.size());
            real[gaps.back()] = real.back();
            gaps.pop_back();
            real.pop_back();
            --diff;
        }

        // total_size -= gaps.size();
        while (gaps.size())
        {
            if (gaps.back() != real.size() - 1)
                real[gaps.back()] = real.back();
            gaps.pop_back();
            real.pop_back();
        }

        PHARE_ASSERT(gaps.size() == 0);
    }
}


template<typename Particles, std::uint8_t impl>
template<auto type>
void PerCellVector<Particles, impl>::sync_gpu_gaps_and_tmp_impl0()
{
    auto const lbox = local_box();

    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync_gpu_gaps_and_tmp::reserve_scan ");
        for (auto const& bix : lbox)
        { // calculate reserve size
            auto& real = particles_(bix);
            if constexpr (any_in(impl_v, 1, 2))
            {
                // push_back is done on device requires over allocation
                // cell_size_(bix) = particles_views_(bix).size();
                // real.resize(particles_views_(bix).size());
                resize(real, particles_views_(bix).size());
            }
            reserve(real, real.size() + add_into_(bix));
        }
    }

    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync_gpu_gaps_and_tmp::add ");
        for (auto const& bix : lbox)
        { // add incoming particles
            auto& real            = particles_(bix);
            auto const& gaps      = gaps_(bix);
            auto const& gaps_size = gap_idx_(bix);
            for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
            {
                auto const& idx     = gaps[gidx];
                auto const& newcell = local_cell(real.iCell(idx));
                // auto const& p = real[idx];

                /*particles_(local_cell(p.iCell())).*/ emplace_back(particles_(newcell), real, idx);
            }
        }
    }

    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync_gpu_gaps_and_tmp::delete ");
        for (auto const& bix : lbox)
        { // delete outgoing particles
            auto const& gaps_size = gap_idx_(bix);
            {
                auto& gaps = gaps_(bix);
                std::sort(gaps.begin(), gaps.begin() + gaps_size, std::greater<>()); // use thrust?
            }
            auto& real       = particles_(bix);
            auto const& gaps = gaps_(bix);
            for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
            {
                auto const& idx = gaps[gidx];
                real.assign(real.size() - 1, idx);
                real.pop_back();
            }
            gap_idx_(bix)  = 0;
            add_into_(bix) = 0;
        }
    }
}

template<typename Particles, std::uint8_t impl>
template<auto type>
void PerCellVector<Particles, impl>::sync_gpu_gaps_and_tmp_impl1()
{
}



template<typename Particles, std::uint8_t impl>
template<auto type>
void PerCellVector<Particles, impl>::sync_gpu_gaps_and_tmp()
{
    // PHARE_LOG_LINE_STR("sync_gpu_gaps_and_tmp " << magic_enum::enum_name(type));
    PHARE_LOG_SCOPE(1, "PerCellVector::sync_gpu_gaps_and_tmp ");

    // if constexpr (impl == 0)
    // {
    sync_gpu_gaps_and_tmp_impl0<type>();
    // }
    // else if constexpr (impl == 1)
    // {
    //     sync_gpu_gaps_and_tmp_impl1<type>();
    // }
    // else
    //     throw std::runtime_error("No impl");
}

template<typename Particles, std::uint8_t impl>
void PerCellVector<Particles, impl>::reset_index_wrapper_map()
{
    resize(p2c_, total_size);

    auto const fill = [](auto p, auto o, auto s, auto b) { std::fill(p + o, p + o + s, *b); };

    std::size_t offset = 0;
    for (auto const& bix : local_box())
    {
        auto const& cs = cell_size_(bix);
        resize(gaps_(bix), cs);
        off_sets_(bix) = offset;

        if (cs)
        {
            if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
            {
                PHARE_WITH_THRUST( //
                    thrust::fill(thrust::device, p2c_.begin() + offset, p2c_.begin() + offset + cs,
                                 *bix));
                PHARE_WITH_THRUST_ELSE(
                    PHARE_LOG_LINE_SS("Thrust not found for PerCellVector<Particles, "
                                      "impl>::reset_index_wrapper_map"); //
                    fill(p2c_.begin(), offset, cs, bix);                 //
                )
            }
            else
                fill(p2c_.begin(), offset, cs, bix);
        }

        offset += cs;
        cap_(bix) = particles_(bix).capacity();
    }
}

template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void PerCellVector<Particles, impl>::sync()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    PHARE_LOG_SCOPE(1, "PerCellVector::sync");
    PHARE_LOG_LINE_STR("sync " << static_cast<std::uint32_t>(PHASE) << " "
                               << magic_enum::enum_name(type));

    auto const lbox = local_box();

    if constexpr (PHASE < 2 and alloc_mode == AllocatorMode::CPU)
        sync_cpu_gaps_and_tmp<type>();

    if constexpr (PHASE < 2 and alloc_mode == AllocatorMode::GPU_UNIFIED)
        sync_gpu_gaps_and_tmp<type>();

    if constexpr (PHASE == 1 || type == ParticleType::Domain)
        trim<ParticleType::Ghost>();

    total_size = 0;
    for (auto const& bix : lbox)
        total_size += (cell_size_(bix) = particles_(bix).size());

    if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync::reset");
        static_assert(impl < 3); // otherwise unhandled
        if constexpr (impl < 2)
        {
            reset_index_wrapper_map();
        }
        else // if (PHASE == 2)
        {
            auto const per_cell = [&](auto& bix) {
                auto const& cs  = cell_size_(bix);
                auto const& cap = particles_(bix).capacity();
                auto& gaps      = gaps_(bix);
                if (gaps.size() < cs)
                {
                    reserve(gaps, cap, false);
                    resize(gaps, cs, false);
                }
                cap_(bix) = cap;
            };

            if constexpr (type == ParticleType::Domain)
                on_domain(per_cell);

            // if constexpr (ParticleType::Domain)
            //     on_ghost_layer(per_cell);
        }
    }

    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync::reset_views");
        reset_views();
    }
}


template<typename Super_>
struct PerCellParticles : public Super_
{
    using Super              = Super_;
    using This               = PerCellParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_cell_particles = typename Super::per_cell_particles;

    auto static constexpr impl_v     = Super::impl_v;
    auto static constexpr alloc_mode = Super::alloc_mode;
    auto static constexpr dimension  = Super::dimension;
    // auto static constexpr layout_mode  = Super::layout_mode;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    // template<std::size_t size>
    // using array_type = AoSParticles<AoSArray<dimension, size>>;

    using Super::local_box;
    using Super::particles_;
    using Super::size;
    // using Super::sync;

    template<typename... Args>
    PerCellParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }



    PerCellParticles(PerCellParticles const& from)            = default;
    PerCellParticles(PerCellParticles&& from)                 = default;
    PerCellParticles& operator=(PerCellParticles&& from)      = default;
    PerCellParticles& operator=(PerCellParticles const& from) = default;

    template<typename T>
    struct iterator_impl;

    template<typename T, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto erase(iterator_impl<T> a, iterator_impl<T> b)
    {
        return Super::erase(a, b);
    }

    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args)
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T&>{*t, args...};
    }
    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, true); }
    auto end() _PHARE_ALL_FN_ { return it(this, true); }

    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }



    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<std::uint32_t, dimension> const& local) const
    {
        std::vector<Particle_t> out;
        out.reserve(sum_from(local, [&](auto const& b) { return particles_(b.toArray()).size(); }));
        for (auto const& b : local)
            std::copy(particles_(b).begin(), particles_(b).end(), std::back_inserter(out));
        return out;
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<int, dimension> const& amr) const
    {
        return (*this)[Super::local_box(amr)];
    }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }


    template<typename _Particle>
    auto& icell_changer(_Particle const& p, std::array<std::uint32_t, dimension> const& cell,
                        std::size_t const& idx,
                        std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            Super::get_vec(Super::gaps_(cell)).emplace_back(idx);
            if (isIn(newcell, Super::ghost_box()))
                Super::get_vec(particles_(Super::local_cell(newcell))).emplace_back(p).iCell()
                    = newcell;
        }
        else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            using Op = Operators<typename Super::SIZE_T, true>;

            // printf("L:%d i %llu ic %u,%u change \n", __LINE__, idx, cell[0], cell[1]);
            Super::gaps_(cell)[Op{Super::gap_idx_(cell)}.increment_return_old()] = idx;

            if (isIn(newcell, Super::ghost_box()))
            {
                auto const nc = Super::local_cell(newcell);
                Op{Super::add_into_(nc)}.increment_return_old();
            }
        }
        else
            throw std::runtime_error("no");

        return *this;
    }


    template<typename T>
    struct index_wrapper;

    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return index_wrapper<This const>{this, s};
    }

    void print() const {}

    void check() const {}

    auto max_size() const
    {
        return max_from(this->particles_,
                        [](auto const& v, auto const& i) { return v.data()[i].size(); });
    }


}; // PerCellParticles<Super>




template<typename OuterSuper>
template<typename T>
struct PerCellParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = Particle<dimension>;
    using pointer           = Particle<dimension>*;
    using reference         = Particle<dimension>&;

    iterator_impl(T& particles_, bool end = false)
        : particles{particles_}
    {
        if (end)
        {
            l = particles.particles_.size() - 1;
            i = particles.size(l);
        }
        else
            inc(); // find first particle
    }
    iterator_impl(iterator_impl&& that)                 = default;
    iterator_impl(iterator_impl const& that)            = default;
    iterator_impl& operator=(iterator_impl&& that)      = default;
    iterator_impl& operator=(iterator_impl const& that) = default;

    auto end() { return iterator_impl{particles, true}; }
    auto& inc()
    {
        i = 0;
        for (; l < particles.particles_.size(); ++l)
            if (*this == this->end() or particles.size(l) > 0)
                break;
        return *this;
    }

    iterator_impl& operator++()
    {
        auto last = end();
        if (*this == last)
            return *this;
        ++i;
        if (*this == last)
            return *this;
        if (i == particles.size(l))
        {
            ++l;
            return inc();
        }
        return *this;
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto operator==(iterator_impl const& that) const { return l == that.l and i == that.i; }
    auto operator!=(iterator_impl const& that) const { return !(*this == that); }
    auto operator<(iterator_impl const& that) const
    {
        // PHARE_LOG_LINE_STR(l << " " << i);
        if (l < that.l)
            return true;
        if (l == that.l)
            return i < that.i;
        return false;
    }

    auto& charge() _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    auto& charge() const _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() const _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() const _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() const _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() const _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    auto& operator*() _PHARE_ALL_FN_ { return particles.data()[l][i]; }
    auto& operator*() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    // auto icell_changer(std::array<int, dimension> const /*newcell*/) { PHARE_LOG_LINE_STR(""); }

    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    T particles;
    std::size_t l = 0, i = 0;
};

template<auto layout_mode, typename Particles>
struct index_wrapper_storage;

#if PHARE_HAVE_THRUST
template<typename Particles>
struct index_wrapper_storage<LayoutMode::SoA, Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using per_cell_particles       = typename std::decay_t<Particles>::per_cell_particles;
    using Particle_t = typename SoAZipParticle_t<per_cell_particles, is_const>::value_type;

    template<typename PerCellParticles_t>
    index_wrapper_storage(PerCellParticles_t p, std::size_t const i) _PHARE_ALL_FN_
        : /*particles{p},*/
          particle{*p, i}
    {
    }

    auto& operator*() _PHARE_ALL_FN_ { return particle; }
    auto& operator*() const _PHARE_ALL_FN_ { return particle; }

    per_cell_particles* particles;
    Particle_t particle;
};
#else

#endif // PHARE_HAVE_THRUST

template<typename Particles>
struct index_wrapper_storage<LayoutMode::AoS, Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using per_cell_particles       = typename std::decay_t<Particles>::per_cell_particles;
    using Particle_t               = typename Particles::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;

    template<typename PerCellParticles_t>
    index_wrapper_storage(PerCellParticles_t p, std::size_t const i) _PHARE_ALL_FN_
        : /*particles{p},*/
          particle{&p->data()[i]}
    {
    }

    auto& operator*() _PHARE_ALL_FN_ { return *particle; }
    auto& operator*() const _PHARE_ALL_FN_ { return *particle; }

    // per_cell_particles* particles;
    Particle_p particle;
};

template<auto layout_mode, typename Particles>
struct index_wrapper_storage
{
    // unused
}; // default;

template<typename T>
struct index_wrapper_super
{
    using per_cell_particles = typename std::decay_t<T>::per_cell_particles;
    using value_type         = index_wrapper_storage<per_cell_particles::layout_mode, T>;
};

template<typename ParticlesSuper>
template<typename T>
struct PerCellParticles<ParticlesSuper>::index_wrapper : public index_wrapper_super<T>::value_type
{
    using outer_t = std::decay_t<T>;
    using Super   = typename index_wrapper_super<T>::value_type;

    auto static constexpr dimension = ParticlesSuper::dimension;
    bool static constexpr is_const  = std::is_const_v<std::remove_reference_t<T>>;
    using Particle_t                = typename outer_t::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;


    index_wrapper(T* pc_particles, std::size_t idx_) _PHARE_ALL_FN_
        : Super{&pc_particles->particles_(cell(pc_particles, idx_)), index(pc_particles, idx_)},
          pc_particles_ptr{pc_particles},
          idx{idx_}
    {
        PHARE_ASSERT((**this).iCell()[0] > -10 and (**this).iCell()[0] < 1000); // bad memory
        if constexpr (dimension > 1)
            PHARE_ASSERT((**this).iCell()[1] > -10 and (**this).iCell()[1] < 1000); // bad memory
        if constexpr (dimension > 2)
            PHARE_ASSERT((**this).iCell()[2] > -10 and (**this).iCell()[2] < 1000); // bad memory
    }

    auto& c() const _PHARE_ALL_FN_ { return cell(pc_particles_ptr, idx); }
    static auto& cell(T* arr, std::size_t const idx) _PHARE_ALL_FN_ { return arr->p2c_[idx]; }
    auto i() const _PHARE_ALL_FN_ { return index(pc_particles_ptr, idx); }
    static auto index(T* arr, std::size_t const idx) _PHARE_ALL_FN_
    {
        return idx - arr->off_sets_(cell(arr, idx));
    }


    auto icell_changer(std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        pc_particles_ptr->icell_changer(**this, c(), i(), newcell);
    }


    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }

    auto& operator*() _PHARE_ALL_FN_ { return *super(); }
    auto& operator*() const _PHARE_ALL_FN_ { return *super(); }

    Particle<dimension> copy() const _PHARE_ALL_FN_
    {
        return {(**this).weight(), (**this).charge(), (**this).iCell(), (**this).delta(),
                (**this).v()};
    }


    T* pc_particles_ptr;
    std::size_t idx = 0;
};




template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles, impl>::sync() _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "PerCellSpan::sync()");

    auto view = *this;
    PHARE_WITH_MKN_GPU({
        auto const& gabox = local_ghost_box_;
        mkn::gpu::GDLauncher{gabox.size()}(
            [=] _PHARE_ALL_FN_() mutable { view.template sync_add_new<PHASE, type>(); });
        mkn::gpu::GDLauncher{gabox.size()}(
            [=] _PHARE_ALL_FN_() mutable { view.template sync_rm_left<PHASE, type>(); });
    })
}

template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type, typename... Args>
void PerCellSpan<Particles, impl>::sync(Args&&... args) _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "PerCellSpan::sync(stream)");

    auto view = *this;
    PHARE_WITH_MKN_GPU({
        auto const& [stream] = std::forward_as_tuple(args...);
        auto const& gabox    = local_ghost_box_;
        mkn::gpu::GDLauncher<true>{gabox.size()}.stream(
            stream, [=] _PHARE_ALL_FN_() mutable { view.template sync_add_new<PHASE, type>(); });
        mkn::gpu::GDLauncher<true>{gabox.size()}.stream(
            stream, [=] _PHARE_ALL_FN_() mutable { view.template sync_rm_left<PHASE, type>(); });
    })
}


template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles, impl>::sync_add_new() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU

    auto const& kidx   = mkn::gpu::idx();
    auto const& bix    = *(local_ghost_box_.begin() + kidx);
    auto const& n_gaps = gap_idx_(bix);
    {
        auto& gaps = gaps_(bix);
        thrust::sort(thrust::seq, gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
    }
    auto& real       = particles_(bix);
    auto& left       = left_(bix);
    auto const& gaps = gaps_(bix);
    for (std::size_t i = 0; i < n_gaps; ++i)
    {
        auto const& gidx = gaps[n_gaps - (1 + i)];
        // auto const& part    = real[gaps[n_gaps - (1 + i)]];
        auto const& newcell = local_cell(real.iCell(gidx));
        auto const& cap     = cap_(newcell);
        auto& nparts        = particles_(newcell);
        auto npidx          = nparts.size();

        while (true)
        {
            if (npidx >= cap)
                return;
            auto inc = npidx + 1;
            auto old = atomicCAS(nparts.size_address(), npidx, inc);
            if (npidx != old)
            {
                ++npidx;
                continue;
            }
            else
                break;
        }

        nparts.assign(real, gidx, npidx);
        // nparts[npidx] = part;
        ++left;
    }
#endif // PHARE_HAVE_MKN_GPU
}


// // working version guaranteed overallocated, no atomics
// template<typename Particles, std::uint8_t impl>
// template<std::uint8_t PHASE, auto type>
// void PerCellSpan<Particles, impl>::sync_add_new() _PHARE_ALL_FN_
// {
// #if PHARE_HAVE_MKN_GPU
//     using Op           = Operators<SIZE_T, true>;
//     auto const& kidx   = mkn::gpu::idx();
//     auto const& bix    = *(local_ghost_box_.begin() + kidx);
//     auto const& n_gaps = gap_idx_(bix);
//     {
//         auto& gaps = gaps_(bix);
//         thrust::sort(thrust::seq, gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
//     }
//     auto& real       = particles_(bix);
//     auto const& gaps = gaps_(bix);
//     auto& left       = left_(bix);
//     for (std::size_t i = 0; i < n_gaps; ++i)
//     {
//         auto const& part   = real[gaps[n_gaps - (1 + i)]];
//         auto const newcell = local_cell(part.iCell());
//         auto& nparts       = particles_(newcell);
//         auto const npidx   = Op{nparts.s}.increment_return_old();
//         PHARE_ASSERT(npidx < cap_(newcell));
//         // if (cap_(newcell) > npidx) // TBC
//         {
//             nparts[npidx] = part;
//             ++left;
//         }
//     }
// #endif // PHARE_HAVE_MKN_GPU
// }


template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles, impl>::sync_rm_left() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU
    auto const& kidx = mkn::gpu::idx();
    auto const& bix  = *(local_ghost_box_.begin() + kidx);
    auto const& gaps = gaps_(bix);
    auto& real       = particles_(bix);
    auto& left       = left_(bix);
    auto& gaps_size  = gap_idx_(bix);
    add_into_(bix) -= left;
    while (left)
    {
        auto const& pidx = gaps[gaps_size - 1];
        real.assign(real.size() - 1, pidx); // real[pidx] = real[real.size() - 1];
        real.pop_back();
        --gaps_size;
        --left;
    }

#endif // PHARE_HAVE_MKN_GPU
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP */
