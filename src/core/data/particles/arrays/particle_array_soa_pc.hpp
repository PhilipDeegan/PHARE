#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SoAPCPRO_HPP

#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SoAPCPRO_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SoAPCPRO_HPP


#include "core/utilities/span.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "particle_array_soa.hpp"

namespace PHARE::core
{

template<std::size_t dim>
class SoAPCSpan
{
    template<typename SoAPCArray>
    auto resolve(SoAPCArray& arr)
    {
        if constexpr (SoAPCArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
        {
            arr.check();
            return *arr.particles_views_;
        }
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = Particle<dim>;

    template<typename SoAPCArray>
    SoAPCSpan(SoAPCArray& arr)
        : particles_{resolve(arr)}
        , size_{arr.size()}
    {
    }

    auto size() const { return size_; }
    void resize(std::size_t s) { particles_.s = s; }
    auto& operator[](std::size_t i) const _PHARE_ALL_FN_ { return particles_[i]; }
    auto& operator[](std::size_t i) _PHARE_ALL_FN_ { return particles_[i]; }

protected:
    NdArrayView<dim, Span<Particle_t>> particles_;
    Span<std::size_t> icell_offsets_;
    std::size_t size_;
};

template<std::size_t dim, auto alloc_mode_ = AllocatorMode::CPU>
class SoAPCVector
{
    using This = SoAPCVector<dim, alloc_mode_>;
    friend class SoAPCSpan<dim>;

public:
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = Particle<dim>;
    using value_type                   = Particle_t;

    template<typename T>
    using vec_helper = std::conditional_t<alloc_mode == AllocatorMode::CPU, //
                                          PHARE::Vector<T, alloc_mode>,
                                          PHARE::BufferedVectorHelper<T, alloc_mode>>;

    using size_t_vec_helper   = vec_helper<std::size_t>;
    using size_t_vector       = typename size_t_vec_helper::vector_t;
    using particle_vec_helper = vec_helper<Particle_t>;
    using particle_vector     = typename particle_vec_helper::vector_t;

    SoAPCVector(box_t const& box = {}, std::size_t ghost_cells = 0)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
    {
        cell_size_.zero();
    }

    SoAPCVector(SoAPCVector const& from)            = default;
    SoAPCVector(SoAPCVector&& from)                 = default;
    SoAPCVector& operator=(SoAPCVector&& from)      = default;
    SoAPCVector& operator=(SoAPCVector const& from) = default;

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
    Particle_t& emplace_back(Particle_t&& p)
    {
        auto& np = particles_(local_cell(p.iCell())).emplace_back(std::forward<Particle_t>(p));
        if constexpr (inc_)
            _inc(np);
        return np;
    }

    void push_back(Particle_t&& p) { emplace_back(std::forward<Particle_t>(p)); }

    void reset() { particles_views_ = generate_spans(); }

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }


    auto& operator()(locell_t const& cell) { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }
    auto local_box() const
    {
        return box_from_zero_to_upper_minus_one(
            ghost_box_.shape().template toArray<std::uint32_t>());
    }
    auto local_box(Box<int, dim> const& from) const
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync(); // after move

    template<auto type>
    void trim();

    void print();

    template<typename T>
    struct index_wrapper;
    // auto operator[](std::size_t const& s) { return index_wrapper<This&>{*this, s}; }
    // auto operator[](std::size_t const& s) const { return index_wrapper<This const&>{*this, s}; }

    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        total_size = 0;
    }


    auto& insert(SoAPCVector const& src);

    auto& insert_domain_from(SoAPCVector const& src);

    auto& icell_changer(Particle_t const& p, std::array<std::uint32_t, dim> const& cell,
                        std::size_t const& idx, std::array<int, dim> const& newcell);

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
        if (this == &that) // just in case
            return;

        gaps_      = that.gaps_;
        cell_size_ = that.cell_size_;
        particles_ = that.particles_;
        box_       = that.box_;
        ghost_box_ = that.ghost_box_;
        // cell_size_ = that.cell_size_;
        // cell_size_ = that.cell_size_;
    }

protected:
    bool static constexpr c_order = true;

    std::size_t ghost_cells_;

    Box<int, dim> box_;
    Box<int, dim> ghost_box_;

    NdArrayVector<dim, particle_vector, c_order, alloc_mode> particles_{local_box().shape()};

    // only used for GPU
    NdArrayVector<dim, Span<Particle_t>, c_order, alloc_mode> particles_views_{local_box().shape()};

    typename vec_helper<std::array<std::uint8_t, dim>>::vector_t p2c_;
    NdArrayVector<dim, std::size_t, c_order, alloc_mode> off_sets_{local_box().shape()};

    NdArrayVector<dim, size_t_vector, c_order, alloc_mode> gaps_{local_box().shape()};

    NdArrayVector<dim, std::size_t, c_order, alloc_mode> cell_size_{local_box().shape()};

    std::size_t total_size = 0;

    auto generate_spans()
    {
        return generate_from<alloc_mode>(
            [&](auto const i) { return make_span(*(particles_.data() + i)); }, particles_);
    }


    template<auto type>
    void sync_gaps_and_tmp();

}; // SoAPCVector<dim, alloc_mode>


template<std::size_t dim, auto alloc_mode>
auto& SoAPCVector<dim, alloc_mode>::icell_changer(Particle_t const& p,
                                                  std::array<std::uint32_t, dim> const& cell,
                                                  std::size_t const& idx,
                                                  std::array<int, dim> const& newcell)
{
    gaps_(cell).emplace_back(idx);
    if (isIn(newcell, ghost_box()))
        particles_(local_cell(newcell)).emplace_back(p).iCell() = newcell;

    return *this;
}


template<std::size_t dim, auto alloc_mode>
template<typename T>
struct SoAPCVector<dim, alloc_mode>::index_wrapper
{
    using outer_t = std::decay_t<T>;

    auto static constexpr dimension = dim;
    bool static constexpr is_const  = std::is_const_v<std::remove_reference_t<T>>;
    using Particle_t                = typename outer_t::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;


    index_wrapper(T particles_, std::size_t idx_)
        : array_{particles_}
        , idx{idx_}
        , p_{&pi(*this)}
    {
        // PHARE_LOG_LINE_STR(idx << " " << p_->iCell()[0] << " "
        //                        << static_cast<std::uint32_t>(Point{c()}[0]));
        // PHARE_ASSERT(p_->iCell()[0] > -100000 and p_->iCell()[0] < 1000000); // bad memory
    }

    auto& c() const { return array_.p2c_[idx]; }
    auto i() const { return idx - array_.off_sets_(c()); }

    template<typename This>
    static auto& pi(This& self)
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

    // auto icell_changer(std::array<int, dimension> const& newcell)
    // {
    //     array_.icell_changer(*p_, c(), i(), newcell);
    // }

    T array_;
    std::size_t idx;
    Particle_p p_;
};


template<std::size_t dim, auto alloc_mode>
auto& SoAPCVector<dim, alloc_mode>::insert(SoAPCVector const& src)
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

template<std::size_t dim, auto alloc_mode>
auto& SoAPCVector<dim, alloc_mode>::insert_domain_from(SoAPCVector const& src)
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

template<std::size_t dim, auto alloc_mode>
template<auto type>
auto& SoAPCVector<dim, alloc_mode>::reserve_ppc(std::size_t const& ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    std::size_t const additional = ppc < 10 ? 2 : ppc * .1; // 10% overallocate

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
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(ppc); });
    }

    return *this;
}

template<std::size_t dim, auto alloc_mode>
void SoAPCVector<dim, alloc_mode>::print()
{
    for (auto const& bix : local_box())
    {
        if (particles_(bix).size() == 1)
        {
            PHARE_LOG_LINE_STR(bix << " " << particles_(bix).size() << " "
                                   << Point{particles_(bix)[0].delta()});
        }
        else
        {
            PHARE_LOG_LINE_STR(bix << " " << particles_(bix).size());
        }
    }
}

template<std::size_t dim, auto alloc_mode>
template<auto type>
void SoAPCVector<dim, alloc_mode>::trim()
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

template<std::size_t dim, auto alloc_mode>
template<auto type>
void SoAPCVector<dim, alloc_mode>::sync_gaps_and_tmp()
{
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



template<std::size_t dim, auto alloc_mode>
template<std::uint8_t PHASE, auto type>
void SoAPCVector<dim, alloc_mode>::sync()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    auto const lbox = local_box();

    if constexpr (PHASE == 1 || type == ParticleType::Domain)
        trim<ParticleType::Ghost>();

    if constexpr (PHASE < 2)
        sync_gaps_and_tmp<type>();

    total_size = 0;
    for (auto const& bix : lbox)
        total_size += (cell_size_(bix) = particles_(bix).size());

    // p2c_.resize(total_size);
    // std::size_t offset = 0;
    // auto amr_box_it    = box_.begin();
    // for (auto const& bix : lbox)
    // {
    //     auto const& cs = cell_size_(bix);
    //     off_sets_(bix) = offset;
    //     if (cs)
    //         std::fill(p2c_.begin() + offset, p2c_.begin() + offset + cs,
    //                   as_local_cell<std::uint8_t>(box_.lower(), (*amr_box_it).toArray()));
    //     offset += cs;
    //     ++amr_box_it;
    // }
    // total_size = offset;

    // PHARE_LOG_LINE_STR("sync " << static_cast<std::uint32_t>(PHASE) << " "
    //                            << static_cast<std::uint32_t>(type) << " " << total_size);
}


template<typename Super_>
struct SoAPCParticles : public Super_
{
    using Super      = Super_;
    using This       = SoAPCParticles<Super>;
    using Particle_t = typename Super::Particle_t;

    auto static constexpr dimension    = Super::dimension;
    auto static constexpr layout_mode  = LayoutMode::AoS;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    template<std::size_t size>
    using array_type = AoSParticles<AoSArray<dimension, size>>;

    using Super::local_box;
    using Super::particles_;
    using Super::size;
    using Super::sync;

    template<typename... Args>
    SoAPCParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }

    SoAPCParticles(SoAPCParticles const& from)            = default;
    SoAPCParticles(SoAPCParticles&& from)                 = default;
    SoAPCParticles& operator=(SoAPCParticles&& from)      = default;
    SoAPCParticles& operator=(SoAPCParticles const& from) = default;

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


    void check() const {}

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

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](std::size_t const& s)
    {
        return (**this)[s];
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](std::size_t const& s) const
    {
        return (**this)[s];
    }

}; // SoAPCParticles<Super>




template<typename OuterSuper>
template<typename T>
struct SoAPCParticles<OuterSuper>::iterator_impl
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

    auto& operator*() _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_STR(l << " " << i);
        return particles.data()[l][i];
    }
    auto& operator*() const _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_STR(l << " " << i);
        return particles.data()[l][i];
    }

    auto icell_changer(std::array<int, dimension> const newcell) { PHARE_LOG_LINE_STR(""); }

    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    T particles;
    std::size_t l = 0, i = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SoAPCPRO_HPP */
