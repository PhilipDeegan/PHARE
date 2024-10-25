#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP


#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/range/range.hpp"


#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace PHARE::core
{

template<std::size_t dim, std::size_t size_>
class AoSArray
{
public:
    auto static constexpr alloc_mode   = AllocatorMode::CPU;
    auto static constexpr storage_mode = StorageMode::ARRAY;
    auto static constexpr dimension    = dim;
    using Particle_t                   = Particle<dim>;
    using container_type               = std::array<Particle_t, size_>;

    auto begin() const { return particles_.begin(); }
    auto begin() { return particles_.begin(); }
    auto end() const { return particles_.end(); }
    auto end() { return particles_.end(); }
    auto constexpr static size() { return size_; }
    auto& operator[](std::size_t i) const { return particles_.data()[i]; }
    auto& operator[](std::size_t i) { return particles_.data()[i]; }

protected:
    container_type particles_;
};



template<std::size_t dim, auto alloc_mode_>
class AoSSpan
{
    static_assert(std::is_same_v<decltype(alloc_mode_), PHARE::AllocatorMode>);

public:
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::SPAN;
    auto static constexpr dimension    = dim;
    using Particle_t                   = Particle<dim>;

    AoSSpan() = default;

    template<typename Container>
    AoSSpan(Container& container)
        : particles_{container.data(), container.size()}
    {
    }

    template<typename Container>
    AoSSpan(Container& container, std::size_t const& beg, std::size_t const& siz)
        : particles_{container.data() + beg, siz}
    {
    }
    template<typename Container>
    AoSSpan(Container& array, std::size_t const& siz)
        : AoSSpan{array, 0, siz}
    {
    }

    auto size() _PHARE_ALL_FN_ { return particles_.size(); }
    void clear() _PHARE_ALL_FN_ { particles_.s = 0; }
    void resize(std::size_t const& s)
    {
        PHARE_ASSERT(s <= particles_.s); // can't be bigger
        particles_.s = s;
    }

    auto& operator[](std::size_t i) const _PHARE_ALL_FN_ { return particles_[i]; }
    auto& operator[](std::size_t i) _PHARE_ALL_FN_ { return particles_[i]; }


    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }

    void pop_back() _PHARE_ALL_FN_ { --particles_.s; }

    auto size_address() _PHARE_ALL_FN_ { return &particles_.s; }


    Span<Particle_t> particles_;
};

template<std::size_t dim, auto alloc_mode_>
class AoSVector
{
    static_assert(std::is_same_v<decltype(alloc_mode_), PHARE::AllocatorMode>);
    using This = AoSVector<dim, alloc_mode_>;

    template<typename Iterator>
    auto check_distance_size_t(Iterator const& start, Iterator const& end)
    {
        auto dist = std::distance(start, end);
        if (dist < 0)
            throw std::runtime_error("Error, number must be postive");
        return static_cast<std::size_t>(dist);
    }


public:
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr dimension    = dim;

    using Particle_t     = Particle<dim>;
    using value_type     = Particle_t;
    using vec_helper     = PHARE::Vector<Particle_t, alloc_mode, 1>;
    using container_type = typename vec_helper::vector_t;

    AoSVector(std::size_t size = 0)
        : particles_{vec_helper::make(size)}
    {
    }

    template<typename Particle_t>
    AoSVector(std::size_t size, Particle_t const& particle)
        : particles_{vec_helper::make(size, particle)}
    {
    }

    template<typename Iterator>
    AoSVector(Iterator start, Iterator end)
        : AoSVector{check_distance_size_t(start, end)}
    {
        std::copy(start, end, particles_.begin());
    }

    AoSVector(AoSVector const& from)            = default;
    AoSVector(AoSVector&& from)                 = default;
    AoSVector& operator=(AoSVector&& from)      = default;
    AoSVector& operator=(AoSVector const& from) = default;



    auto size() const { return particles_.size(); }
    auto capacity() const { return particles_.capacity(); }
    void clear() { particles_.clear(); }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }
    auto& operator[](std::size_t i) const _PHARE_ALL_FN_ { return particles_.data()[i]; }
    auto& operator[](std::size_t i) _PHARE_ALL_FN_ { return particles_.data()[i]; }
    bool operator==(This const& that) const { return (this->particles_ == that.particles_); }


    void pop_back() { particles_.pop_back(); }

    auto back() { return particles_.back(); }
    auto front() { return particles_.front(); }

    template<typename Iterator>
    auto erase(Iterator first, Iterator last)
    {
        // should we erase particles indexes associated with these iterators from the cellmap?
        // probably it does not matter if not. The reason is that
        // particles erased from the particlearray are so because they left
        // the patch cells to an outside cell.
        // But in principle that cell will never be accessed because it is outside the patch.
        // The only thing "bad" if these indexes are not deleted is that the
        // size of the cellmap becomes unequal to the size of the particleArray.
        // but  ¯\_(ツ)_/¯
        return particles_.erase(particles_.begin() + first.curr_pos,
                                particles_.begin() + last.curr_pos);
    }

    Particle_t& emplace_back() { return get_vec(particles_).emplace_back(); }

    Particle_t& emplace_back(Particle_t&& p)
    {
        return get_vec(particles_).emplace_back(std::forward<Particle_t>(p));
    }

    Particle_t& emplace_back(Particle_t const& p) { return get_vec(particles_).emplace_back(p); }

    void emplace_back(This const& src, std::size_t const& idx) { emplace_back(src[idx]); }
    void emplace_back(This const& src)
    {
        for (auto const& p : src.particles_)
            emplace_back(p);
    }

    template<typename... Args>
    Particle_t& emplace_back(double const& weight, Args const&... args)
    {
        return get_vec(particles_).emplace_back(weight, args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(double const& weight, Args&&... args)
    {
        return get_vec(particles_).emplace_back(weight, args...);
    }


    void push_back(Particle_t const& p) { get_vec(particles_).push_back(p); }
    void push_back(Particle_t&& p) { get_vec(particles_).push_back(std::forward<Particle_t>(p)); }

    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }

protected:
    container_type particles_;

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
};


template<typename Super_>
struct AoSParticles : public Super_
{
    using Super      = Super_;
    using This       = AoSParticles<Super>;
    using Particle_t = typename Super::Particle_t;

    auto static constexpr dimension    = Super::dimension;
    auto static constexpr layout_mode  = LayoutMode::AoS;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Span_t = AoSParticles<AoSSpan<dimension, alloc_mode>>;
    friend struct AoSParticles<AoSSpan<dimension, alloc_mode>>;


    template<std::size_t size>
    using array_type = AoSParticles<AoSArray<dimension, size>>;

    using Super::particles_;

    template<typename... Args>
    AoSParticles(Args&&... args) _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    AoSParticles(AoSParticles const& from)            = default;
    AoSParticles(AoSParticles&& from)                 = default;
    AoSParticles& operator=(AoSParticles&& from)      = default;
    AoSParticles& operator=(AoSParticles const& from) = default;


    template<typename T>
    struct iterator_impl;

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto erase(iterator_impl<This*> a, iterator_impl<This*> b)
    {
        return Super::erase(a, b);
    }

    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args) _PHARE_ALL_FN_
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T*>{t, args...};
    }
    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, size()); }
    auto end() _PHARE_ALL_FN_ { return it(this, size()); }


    auto& weight(std::size_t i) const { return particles_[i].weight(); }
    auto& weight(std::size_t i) { return particles_[i].weight(); }

    auto& charge(std::size_t i) const { return particles_[i].charge(); }
    auto& charge(std::size_t i) _PHARE_ALL_FN_ { return particles_[i].charge(); }

    auto& iCell(std::size_t i) const _PHARE_ALL_FN_ { return particles_[i].iCell(); }
    auto& iCell(std::size_t i) _PHARE_ALL_FN_ { return particles_[i].iCell(); }

    auto& delta(std::size_t i) const _PHARE_ALL_FN_ { return particles_[i].delta(); }
    auto& delta(std::size_t i) _PHARE_ALL_FN_ { return particles_[i].delta(); }

    auto& v(std::size_t i) const { return particles_[i].v(); }
    auto& v(std::size_t i) _PHARE_ALL_FN_ { return particles_[i].v(); }


    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }


    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto capacity() const
    {
        return particles_.capacity();
    }
    auto size() const _PHARE_ALL_FN_ { return particles_.size(); }

    template<typename IndexRange, typename Predicate>
    auto partition(IndexRange&& range, Predicate&& pred)
    {
        return std::partition(range.begin(), range.end(), pred);
    }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto& vector()
    {
        return particles_;
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto& vector() const
    {
        return particles_;
    }

    // does not swap cellmap (never did)
    void swap(This& that) { std::swap(this->particles_, that.particles_); }
    void swap(std::size_t const& a, std::size_t const& b)
    {
        if (a == b)
            return;

        std::swap(particles_[a], particles_[b]);
    }

    void replace_from(This const& that)
    {
        if (this == &that) // just in case
            return;
        this->resize(that.size());
        std::copy(that.begin(), that.end(), this->begin());
    }

    void check() const {}

    template<typename _Particles>
    void assign(_Particles const& src, std::size_t const& idx,
                std::size_t const& dst) _PHARE_ALL_FN_
    {
        particles_[dst] = src[idx];
    }

    void assign(Particle_t const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        particles_[dst] = src;
    }

    void assign(std::size_t const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        return assign(particles_[src], dst);
    }
};



template<typename Super_>
class AoSMappedParticles : public Super_
{
    using Super = Super_;
    using This  = AoSMappedParticles<Super>;
    using Super::particles_;

    template<typename _Super>
    friend class AoSMappedParticles;

public:
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;

    using box_t      = Box<int, dimension>;
    using CellMap_t  = CellMap<dimension, int>;
    using Particle_t = typename Super::Particle_t;

    template<std::size_t size>
    using array_type = typename Super::template array_type<size>;


    AoSMappedParticles(std::size_t size)
        : Super_(size)
        , box_{}
        , cellMap_{box_}
    {
    }

    template<typename Particle_t>
    AoSMappedParticles(std::size_t size, Particle_t const& from)
        : Super_(size, from)
        , box_{}
        , cellMap_{box_}
    {
    }

    AoSMappedParticles(box_t box = {}, std::size_t size = 0)
        : Super_(size)
        , box_{box}
        , cellMap_{box_}
    {
    }

    template<typename Particle_t>
    AoSMappedParticles(box_t box, std::size_t size, Particle_t const& from)
        : Super_(size, from)
        , box_{box}
        , cellMap_{box_}
    {
        PHARE_ASSERT(box_.size() > 0);
    }


    AoSMappedParticles(AoSMappedParticles&&)                 = default;
    AoSMappedParticles(AoSMappedParticles const&)            = default;
    AoSMappedParticles& operator=(AoSMappedParticles&&)      = default;
    AoSMappedParticles& operator=(AoSMappedParticles const&) = default;

    template<typename ParticleArray_t>
    AoSMappedParticles(ParticleArray_t& array)
        : Super{array}
        , box_{array.box()}
        , cellMap_{array.cellMap_}
    {
    }


    void clear()
    {
        Super::clear();
        cellMap_.clear();
    }

    template<typename IndexRange>
    auto erase(IndexRange& range)
    {
        cellMap_.erase(particles_, range);
    }
    template<typename IndexRange>
    auto erase(IndexRange&& range)
    {
        // TODO move ctor for range?
        cellMap_.erase(std::forward<IndexRange>(range));
    }

    template<typename It>
    auto erase(It a, It b)
    {
        return Super::erase(a, b);
    }


    auto& emplace_back()
    {
        auto& part = Super::emplace_back();
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }


    auto& emplace_back(Particle_t&& p)
    {
        auto& part = Super::emplace_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    auto& emplace_back(Particle_t const& p)
    {
        auto& part = Super::emplace_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    template<typename... Args>
    auto& emplace_back(Args const&... args)
    {
        auto& part = Super::emplace_back(args...);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }
    template<typename... Args>
    auto& emplace_back(Args&&... args)
    {
        auto& part = Super::emplace_back(args...);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    void push_back(Particle_t const& p)
    {
        Super::push_back(p);
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void push_back(Particle_t&& p)
    {
        Super::push_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
    }


    void map_particles() const { cellMap_.add(particles_); }
    void empty_map() { cellMap_.empty(); }

    template<typename... Args>
    auto nbr_particles_in(Args&&... args) const
    {
        return cellMap_.size(args...);
    }

    void export_particles(box_t const& box, This& dest) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles");
        cellMap_.export_to(box, *this, dest);
    }

    template<typename Dest, typename Fn>
    void export_particles(box_t const& box, Dest& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn)");
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (box, vector, Fn)");
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(This& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn,vector)");
        cellMap_.export_if(*this, dest, std::forward<Predicate>(pred));
    }


    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t particleIndex)
    {
        auto oldCell                      = particles_[particleIndex].iCell();
        particles_[particleIndex].iCell() = newCell;
        if (!box_.isEmpty())
        {
            cellMap_.update(particles_, particleIndex, oldCell);
        }
    }


    template<typename IndexRange, typename Predicate>
    auto partition(IndexRange&& range, Predicate&& pred)
    {
        return cellMap_.partition(range, std::forward<Predicate>(pred));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cellMap_.print(cell);
    }


    void sortMapping() const { cellMap_.sort(); }

    bool is_mapped() const
    {
        if (particles_.size() != cellMap_.size())
        {
            throw std::runtime_error("particle array not mapped, map.size() != array.size()");
        }
        for (std::size_t pidx = 0; pidx < particles_.size(); ++pidx)
        {
            auto const& p = particles_[pidx];
            auto& icell   = p.iCell;
            auto l        = cellMap_.list_at(icell);
            if (!l)
                throw std::runtime_error("particle cell not mapped");
            auto& ll = l->get();
            if (!ll.is_indexed(pidx))
                throw std::runtime_error("particle not indexed");
        }
        return true;
    }


    void check() const
    {
        if (particles_.size() != cellMap_.size())
        {
            PHARE_LOG_LINE_STR(particles_.size());
            PHARE_LOG_LINE_STR(cellMap_.size());
        }

        core::abort_if(particles_.size() != cellMap_.size());
        // throw std::runtime_error("particle array not mapped, map.size() != array.size()");
    }

    auto& box() const { return box_; }


    void replace_from(This const& that)
    {
        if (this == &that) // just in case
            return;
        Super::replace_from(that);
        this->box_     = that.box_;
        this->cellMap_ = that.cellMap_;
    }


protected:
    box_t box_{};
    mutable CellMap_t cellMap_;
};




template<typename OuterSuper>
template<typename T>
struct AoSParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = Particle<dimension>;
    using pointer           = Particle<dimension>*;
    using reference         = Particle<dimension>&;


    iterator_impl(T& particles_, std::size_t const s = 0) _PHARE_ALL_FN_ : particles{particles_},
                                                                           curr_pos{s}
    {
    }
    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    iterator_impl& operator=(iterator_impl&& that)      = default;
    iterator_impl& operator=(iterator_impl const& that) = default;


    auto& operator++()
    {
        ++curr_pos;
        return *this;
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto& operator+=(std::int64_t i) _PHARE_ALL_FN_
    {
        curr_pos += i;
        return *this;
    }


    auto& operator--() _PHARE_ALL_FN_
    {
        --curr_pos;
        return *this;
    }
    auto operator+(std::int64_t i) const _PHARE_ALL_FN_
    {
        auto copy = *this;
        copy.curr_pos += i;
        return copy;
    }
    auto operator-(std::int64_t i) const
    {
        auto copy = *this;
        copy.curr_pos -= i;
        return copy;
    }


    auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
    auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
    auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
    auto operator-(iterator_impl const& that) const { return curr_pos - that.curr_pos; }

    auto& weight() _PHARE_ALL_FN_ { return deref(particles)[curr_pos].weight(); }
    auto& weight() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos].weight(); }

    auto& charge() _PHARE_ALL_FN_ { return deref(particles)[curr_pos].charge(); }
    auto& charge() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos].charge(); }

    auto& iCell() _PHARE_ALL_FN_ { return deref(particles)[curr_pos].iCell(); }
    auto& iCell() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos].iCell(); }

    auto& delta() _PHARE_ALL_FN_ { return deref(particles)[curr_pos].delta(); }
    auto& delta() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos].delta(); }

    auto& v() _PHARE_ALL_FN_ { return deref(particles)[curr_pos].v(); }
    auto& v() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos].v(); }

    auto& operator*() _PHARE_ALL_FN_
    {
        // PHARE_ASSERT(curr_pos < particles->size());
        return deref(particles)[curr_pos];
    }
    auto& operator*() const _PHARE_ALL_FN_
    {
        // PHARE_ASSERT(curr_pos < particles->size());
        return deref(particles)[curr_pos];
    }

    auto& operator[](std::size_t const i) _PHARE_ALL_FN_ { return deref(particles)[curr_pos + i]; }
    auto& operator[](std::size_t const i) const _PHARE_ALL_FN_
    {
        return deref(particles)[curr_pos + i];
    }

    auto copy() const _PHARE_ALL_FN_ { return deref(particles)[curr_pos]; }

    T particles;
    std::size_t curr_pos;
};


} // namespace PHARE::core




#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP */
