#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AoSPCPRO_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AoSPCPRO_HPP

#include "core/operators.hpp"
#include "core/utilities/span.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "particle_array_aos.hpp"

#include <unordered_set>

namespace PHARE::core
{

template<std::size_t dim, auto alloc_mode_, std::uint8_t impl_ = 0>
class AoSPCSpan
{
protected:
    using SIZE_T = unsigned long long int;

public:
    auto static constexpr impl_v     = impl_;
    auto static constexpr alloc_mode = alloc_mode_;
    using This                       = AoSPCSpan<dim, alloc_mode, impl_v>;

private:
    using lobox_t  = Box<std::uint32_t, dim>;
    using locell_t = std::array<std::uint32_t, dim>;

    template<typename AoSPCArray>
    auto resolve(AoSPCArray& arr)
    {
        if constexpr (AoSPCArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
        {
            arr.check();
            return *arr.particles_views_;
        }
    }

    template<typename AoSPCArray>
    auto resolve_gaps(AoSPCArray& arr)
    {
        if constexpr (AoSPCArray::storage_mode == StorageMode::SPAN)
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

    template<typename AoSPCArray>
    AoSPCSpan(AoSPCArray& arr)
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , off_sets_{arr.off_sets_}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cap_{arr.cap_}
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

    void resize(std::size_t s) { particles_.s = s; }

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
    auto local_box() const _PHARE_ALL_FN_
    {
        return box_from_zero_to_upper_minus_one(
            safe_box_.shape().template toArray<std::uint32_t>());
    }
    auto local_box(Box<int, dim> const& from) const _PHARE_ALL_FN_
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync() _PHARE_ALL_FN_;

    template<typename Vec>
    void cap(Vec& v)
    {
        if constexpr (impl_v == 1)
        {
            for (auto& p : v())
                p.resize(p.capacity());
            _cap = true;
        }
    }

protected:
    template<auto type>
    void sync_gpu_gaps_and_tmp_impl0() _PHARE_ALL_FN_;

    template<auto type>
    void sync_gpu_gaps_and_tmp_impl1() _PHARE_ALL_FN_;

    NdArrayView<dim, Span<Particle_t>> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> off_sets_, gap_idx_, add_into_, cap_;
    Span<std::array<std::uint32_t, dim>> p2c_;
    std::size_t size_;

    Box<int, dim> box_, ghost_box_, safe_box_;
    lobox_t local_ghost_box_;

    bool _cap = false;
};


template<std::size_t dim, auto alloc_mode_, std::uint8_t impl_ = 0>
class AoSPCVector
{
    using This = AoSPCVector<dim, alloc_mode_, impl_>;
    friend class AoSPCSpan<dim, alloc_mode_, impl_>;

protected:
    using SIZE_T = unsigned long long int; // cuda issues

public:
    auto static constexpr impl_v       = impl_;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = Particle<dim>;
    using value_type                   = Particle_t;

    // template<typename T>
    // using vec_helper = std::conditional_t<alloc_mode == AllocatorMode::CPU, //
    //                                       PHARE::Vector<T, alloc_mode>,
    //                                       PHARE::BufferedVectorHelper<T, alloc_mode>>;
    template<typename T>
    using vec_helper = PHARE::Vector<T, alloc_mode>;

    using size_t_vec_helper   = vec_helper<std::size_t>;
    using size_t_vector       = typename size_t_vec_helper::vector_t;
    using particle_vec_helper = vec_helper<Particle_t>;
    using particle_vector     = typename particle_vec_helper::vector_t;

    AoSPCVector(box_t const& box = {}, std::size_t ghost_cells = 0)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
        , safe_box_{grow(box, ghost_cells + 1)}
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        cap_.zero();
    }

    AoSPCVector(AoSPCVector const& from)            = default;
    AoSPCVector(AoSPCVector&& from)                 = default;
    AoSPCVector& operator=(AoSPCVector&& from)      = default;
    AoSPCVector& operator=(AoSPCVector const& from) = default;

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

    void _inc(Particle_t& p)
    {
        ++total_size;
        ++cell_size_(local_cell(p.iCell()));
    }

    template<bool inc_ = true>
    Particle_t& emplace_back(Particle_t const& p)
    {
        auto& np = particles_(local_cell(p.iCell())).emplace_back(p);
        if constexpr (inc_)
            _inc(np);
        return np;
    }

    void push_back(Particle_t&& p) { emplace_back(p); }
    void push_back(Particle_t const& p) { emplace_back(p); }

    void reset_views()
    {
        particles_views_ = generate_from<alloc_mode>(
            [&](auto const i) { return make_span(*(particles_.data() + i)); }, particles_);
        gap_views_ = generate_from<alloc_mode>(
            [&](auto const i) { return make_span(*(gaps_.data() + i)); }, gaps_);
    }

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }


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


    auto& insert(AoSPCVector const& src);

    auto& insert_domain_from(AoSPCVector const& src);


    void check() const
    {
        for (auto const& particles : particles_)
            for (auto const& particle : particles)
            {
                PHARE_ASSERT(particle.iCell()[0] < 1000);
            }
    }

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

        // cell_size_ = that.cell_size_;
        // cell_size_ = that.cell_size_;
    }

protected:
    bool static constexpr c_order = true;

    std::size_t ghost_cells_;

    Box<int, dim> box_, ghost_box_, safe_box_;

    NdArrayVector<dim, particle_vector, c_order, alloc_mode> particles_{local_box().shape()};

    // only used for GPU
    NdArrayVector<dim, Span<Particle_t>, c_order, alloc_mode> particles_views_{local_box().shape()};

    typename vec_helper<std::array<std::uint32_t, dim>>::vector_t p2c_;
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> off_sets_{local_box().shape()};

    NdArrayVector<dim, size_t_vector, c_order, alloc_mode> gaps_{local_box().shape()};
    NdArrayVector<dim, Span<std::size_t>, c_order, alloc_mode> gap_views_{local_box().shape()};

    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> gap_idx_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> add_into_{local_box().shape()};
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

}; // AoSPCVector<dim, alloc_mode>



template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
auto& AoSPCVector<dim, alloc_mode, impl>::insert(AoSPCVector const& src)
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

template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
auto& AoSPCVector<dim, alloc_mode, impl>::insert_domain_from(AoSPCVector const& src)
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

    on_box_list(local_box(box()).remove(shrink(local_box(box()), 1)), [&](auto const& bix) {
        auto& from = src(bix);
        auto& to   = (*this)(bix);
        added += from.size();
        to.reserve(to.size() + from.size());
        std::copy(from.begin(), from.end(), std::back_inserter(to));
    });

    if (added)
        sync<2>();

    return *this;
}

template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
auto& AoSPCVector<dim, alloc_mode, impl>::reserve_ppc(std::size_t const& ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    std::size_t const additional = ppc < 10 ? 5 : ppc * .3; // 30% overallocate

    auto const on_box = [&](auto&& box, auto&& fn) {
        for (auto const& bix : box)
            fn(bix);
    };

    auto const on_box_list = [&](auto&& boxlist, auto&& fn) {
        for (auto const& box : boxlist)
            on_box(box, fn);
    };

    auto const on_domain    = [&](auto&& fn) { on_box(local_box(box()), fn); };
    auto const on_ghost_box = [&](auto&& fn) { on_box(local_box(), fn); };
    auto const on_ghost_layer
        = [&](auto&& fn) { on_box_list(local_box().remove(local_box(box())), fn); };
    auto const on_ghost_layer_plus_1_domain
        = [&](auto&& fn) { on_box_list(local_box().remove(shrink(local_box(box()), 1)), fn); };

    if constexpr (type == ParticleType::Domain)
    {
        on_domain([&](auto const& bix) { particles_(bix).reserve(ppc + additional); });
        on_ghost_box([&](auto const& bix) { gaps_(bix).reserve(additional); });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(additional); });
    }

    if constexpr (type == ParticleType::Ghost)
    {
        on_ghost_layer_plus_1_domain([&](auto const& bix) {
            particles_(bix).reserve(additional);
            gaps_(bix).reserve(additional);
        });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(ppc + additional); });
    }

    return *this;
}



template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCVector<dim, alloc_mode, impl>::trim()
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

template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCVector<dim, alloc_mode, impl>::sync_cpu_gaps_and_tmp()
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


template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCVector<dim, alloc_mode, impl>::sync_gpu_gaps_and_tmp_impl0()
{
    auto const lbox = local_box();

    {
        PHARE_LOG_SCOPE(1, "AoSPCVector::sync_gpu_gaps_and_tmp::reserve_scan ");
        for (auto const& bix : lbox)
        { // calculate reserve size
            auto& real = particles_(bix);
            if constexpr (impl_v == 1)
            {
                // push_back is done on device requires over allocation
                real.resize(particles_views_(bix).size());
            }
            real.reserve(real.size() + add_into_(bix));
        }
    }

    {
        PHARE_LOG_SCOPE(1, "AoSPCVector::sync_gpu_gaps_and_tmp::add ");
        for (auto const& bix : lbox)
        { // add incoming particles
            auto& real            = particles_(bix);
            auto const& gaps      = gaps_(bix);
            auto const& gaps_size = gap_idx_(bix);
            for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
            {
                auto const& idx = gaps[gidx];
                PHARE_ASSERT(idx >= 0 and idx < 1000);
                auto const& p = real[idx];
                PHARE_ASSERT(Point{p.iCell()} != bix);
                particles_(local_cell(p.iCell())).emplace_back(p);
            }
        }
    }

    {
        PHARE_LOG_SCOPE(1, "AoSPCVector::sync_gpu_gaps_and_tmp::delete ");
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
                PHARE_ASSERT(idx >= 0 and idx < 1000);
                if (idx < real.size() - 1)
                    real[idx] = real.back();
                real.pop_back();
            }
            gap_idx_(bix)  = 0;
            add_into_(bix) = 0;
        }
    }
}

template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCVector<dim, alloc_mode, impl>::sync_gpu_gaps_and_tmp_impl1()
{
    // nothing
}



template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCVector<dim, alloc_mode, impl>::sync_gpu_gaps_and_tmp()
{
    // PHARE_LOG_LINE_STR("sync_gpu_gaps_and_tmp " << magic_enum::enum_name(type));
    PHARE_LOG_SCOPE(1, "AoSPCVector::sync_gpu_gaps_and_tmp ");

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

template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void AoSPCVector<dim, alloc_mode, impl>::sync()
{
    PHARE_LOG_LINE_STR("sync " << static_cast<std::uint32_t>(PHASE) << " " << magic_enum::enum_name(type));

    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

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
        p2c_.resize(total_size);
        std::size_t offset = 0;
        for (auto const& bix : lbox)
        {
            auto const& cs = cell_size_(bix);
            gaps_(bix).resize(cs, 0);
            off_sets_(bix) = offset;
            if (cs)
                std::fill(p2c_.begin() + offset, p2c_.begin() + offset + cs, *bix);
            offset += cs;

            cap_(bix) = particles_(bix).capacity();
        }
    }

    reset_views();

    if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        for (auto const& p2c : p2c_)
        {
            PHARE_ASSERT(particles_(p2c).size());
        }
    }
}


template<typename Super_>
struct AoSPCParticles : public Super_
{
    using Super      = Super_;
    using This       = AoSPCParticles<Super>;
    using Particle_t = typename Super::Particle_t;

    auto static constexpr impl_v       = Super::impl_v;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr layout_mode  = LayoutMode::AoS;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    template<std::size_t size>
    using array_type = AoSParticles<AoSArray<dimension, size>>;

    // using Super::local_box;
    using Super::particles_;
    using Super::size;
    // using Super::sync;

    template<typename... Args>
    AoSPCParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }

    AoSPCParticles(AoSPCParticles const& from)            = default;
    AoSPCParticles(AoSPCParticles&& from)                 = default;
    AoSPCParticles& operator=(AoSPCParticles&& from)      = default;
    AoSPCParticles& operator=(AoSPCParticles const& from) = default;

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


    auto& icell_changer(Particle_t const& p, std::array<std::uint32_t, dimension> const& cell,
                        std::size_t const& idx,
                        std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            Super::gaps_(cell).emplace_back(idx);
            if (isIn(newcell, Super::ghost_box()))
                particles_(Super::local_cell(newcell)).emplace_back(p).iCell() = newcell;
        }
        else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            using Op = Operators<typename Super::SIZE_T, true>;
            
            auto const kidx  = mkn::gpu::idx();
            auto const nc = Super::local_cell(newcell);
            auto const oc = Super::local_cell(p.iCell());
            __threadfence();
            auto const nidx = Op{Super::gap_idx_(cell)}.increment_return_old();
            // printf("L:%d k %u c %u,%u oc %u,%u nc %u,%u i %lu x %lu ni %llu\n", __LINE__, kidx, cell[0], cell[1], oc[0], oc[1], nc[0], nc[1], p.id, idx, nidx);

            __threadfence();
            // __syncthreads();
            Op{Super::add_into_(nc)}.increment_return_old();
            // __syncthreads();
            Super::gaps_(cell)[nidx] = idx;
            // __syncthreads();

            // Super::gaps_(cell)[idx] = 1;
        }
        else
            throw std::runtime_error("no");

        return *this;
    }


    template<typename T>
    struct index_wrapper;
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return index_wrapper<This&>{*this, s}; }
    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return index_wrapper<This const&>{*this, s};
    }

    void print() const {}

    void check() const
    {
        if (size() == 0)
            return;

        if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            auto& c = Super::p2c_[0];
            PHARE_ASSERT(&*(*this)[0] == &*this->begin());
        }
    }


}; // AoSPCParticles<Super>




template<typename OuterSuper>
template<typename T>
struct AoSPCParticles<OuterSuper>::iterator_impl
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

    auto icell_changer(std::array<int, dimension> const newcell) { PHARE_LOG_LINE_STR(""); }

    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    T particles;
    std::size_t l = 0, i = 0;
};



template<typename Super_>
template<typename T>
struct AoSPCParticles<Super_>::index_wrapper
{
    using outer_t = std::decay_t<T>;

    auto static constexpr dimension = Super_::dimension;
    bool static constexpr is_const  = std::is_const_v<std::remove_reference_t<T>>;
    using Particle_t                = typename outer_t::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;


    index_wrapper(T particles_, std::size_t idx_) _PHARE_ALL_FN_ : array_{particles_},
                                                                   idx{idx_},
                                                                   p_{&pi(*this)}
    {
        PHARE_ASSERT(p_);
        PHARE_ASSERT(p_->iCell()[0] > -10 and p_->iCell()[0] < 1000); // bad memory
        if constexpr (dimension > 1)
        {
            PHARE_ASSERT(p_->iCell()[1] > -10 and p_->iCell()[1] < 1000); // bad memory
        }
    }

    auto& c() const _PHARE_ALL_FN_ { return array_.p2c_[idx]; }
    auto i() const _PHARE_ALL_FN_ { return idx - array_.off_sets_(c()); }

    template<typename This>
    static auto& pi(This& self) _PHARE_ALL_FN_
    {
        return self.array_.particles_(self.c())[self.i()];
    }

    auto& charge() _PHARE_ALL_FN_ { return p_->charge(); }
    auto& delta() _PHARE_ALL_FN_ { return p_->delta(); }
    auto& iCell() _PHARE_ALL_FN_ { return p_->iCell(); }
    auto& weight() _PHARE_ALL_FN_ { return p_->weight(); }
    auto& v() _PHARE_ALL_FN_ { return p_->v(); }

    auto& charge() const _PHARE_ALL_FN_ { return p_->charge(); }
    auto& delta() const _PHARE_ALL_FN_ { return p_->delta(); }
    auto& iCell() const _PHARE_ALL_FN_ { return p_->iCell(); }
    auto& weight() const _PHARE_ALL_FN_ { return p_->weight(); }
    auto& v() const _PHARE_ALL_FN_ { return p_->v(); }

    auto& operator*() { return *p_; }
    auto& operator*() const { return *p_; }

    auto icell_changer(std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        array_.icell_changer(*p_, c(), i(), newcell);
    }

    auto copy() const { return *p_; }

    T array_;
    std::size_t idx = 0;
    Particle_p p_   = nullptr;
};


template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCSpan<dim, alloc_mode, impl>::sync_gpu_gaps_and_tmp_impl0() _PHARE_ALL_FN_
{
    // nothing
}



template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<auto type>
void AoSPCSpan<dim, alloc_mode, impl>::sync_gpu_gaps_and_tmp_impl1() _PHARE_ALL_FN_
{
    
        using Op = Operators<SIZE_T, true>;

        __syncthreads();
    
        auto const& kidx  = mkn::gpu::idx();
        auto const& gabox = local_ghost_box_;
        // printf("L:%d k %u n %lu \n", __LINE__, kidx, gabox.size());
        // __threadfence();

        std::size_t pop = 0;
        // printf("k %u n %u \n", kidx, n_gaps);

        // Op{Super::add_into_(Super::local_cell(newcell))}.increment_return_old();
        // Super::gaps_(cell)[Op{Super::gap_idx_(cell)}.increment_return_old()] = idx;

        auto const copy_out = [&]() {
            auto const& bix    = *(gabox.begin() + kidx);
            auto const& n_gaps = gap_idx_(bix);

            // printf("L:%d k %u n %llu bix %u,%u \n", __LINE__, kidx, n_gaps, bix[0], bix[1]);
            // __threadfence();
            
            if (n_gaps == 0)
                return;
            
            std::size_t const it = 1;
            {
                auto& gaps = gaps_(bix);
                thrust::sort(thrust::seq, gaps.data(), gaps.data() + n_gaps/*, std::greater<>()*/);
            }
            auto& real       = particles_(bix);
            auto const& gaps = gaps_(bix);
            for (std::size_t i = 0; i < n_gaps; ++i)
            {
                auto const& part   = real[gaps[n_gaps - (1+i)]];
                auto const newcell = local_cell(part.iCell());
                if constexpr (dim > 1)
                {
                    // printf("L:%d k %u n %lu bix %u,%u newcell %u,%u i %lu \n", __LINE__, kidx, gaps[i], bix[0], bix[1], newcell[0], newcell[1], part.id);
                }
                auto& nparts = particles_(newcell);
                __threadfence();
                auto const npidx = Op{nparts.s}.increment_return_old();
                // if (cap_(newcell) > npidx)
                {
                    nparts[npidx] = part;
                    ++pop;
                }
            }
        };

        if (kidx < gabox.size())
            copy_out();

        // __threadfence();
        __syncthreads();

        auto const sync_left = [&]() {            
            auto const& bix  = *(gabox.begin() + kidx);
            
            // printf("L:%d k %u p %lu bix %u,%u \n", __LINE__, kidx, pop, bix[0], bix[1]);
            
            auto& real       = particles_(bix);
            auto const& gaps = gaps_(bix);
            auto& gaps_size  = gap_idx_(bix);
            for (std::size_t i = 0; i < pop; ++i)
            {
                auto const& pidx = gaps[gaps_size - 1];
                // if (pidx != real.s - 1)
                real[pidx] = real[real.s - 1];
                --real.s;
                --gaps_size;
            }
        };

        if (kidx < gabox.size() && pop)
            sync_left();

        // __syncthreads();
    
}



template<std::size_t dim, auto alloc_mode, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void AoSPCSpan<dim, alloc_mode, impl>::sync() _PHARE_ALL_FN_
{
    if constexpr (impl == 0)
        sync_gpu_gaps_and_tmp_impl0<type>();
    else if constexpr (impl == 1)
        sync_gpu_gaps_and_tmp_impl1<type>();
}




} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AoSPCPRO_HPP */
