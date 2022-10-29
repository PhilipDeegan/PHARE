#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP


#include <vector>
#include <cstddef>
#include <utility>

#include "particle.hpp"

#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/utilities/iterators.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/range/range.hpp"

namespace PHARE::core
{
template<std::size_t dim, std::size_t size_>
class AoSArray
{
public:
    static constexpr auto dimension = dim;
    using Particle_t                = Particle<dim>;
    using container_type            = std::array<Particle_t, size_>;


    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

    auto begin() const { return particles_.begin(); }
    auto begin() { return particles_.begin(); }

    auto end() const { return particles_.end(); }
    auto end() { return particles_.end(); }

    auto constexpr static size() { return size_; }

protected:
    container_type particles_;
};

template<std::size_t dim>
class AoSVector
{
    template<typename Iterator>
    auto check_distance_size_t(Iterator const& start, Iterator const& end)
    {
        auto dist = std::distance(start, end);
        if (dist < 0)
            throw std::runtime_error("Error, number must be postive");
        return static_cast<std::size_t>(dist);
    }

    using This = AoSVector<dim>;

public:
    static constexpr auto dimension = dim;

    using Particle_t     = Particle<dim>;
    using value_type     = Particle_t;
    using container_type = std::vector<Particle_t>;


    AoSVector(std::size_t size = 0)
        : particles_(size)
    {
    }

    template<typename Particle_t>
    AoSVector(std::size_t size, Particle_t&& particle)
        : particles_(size, particle)
    {
    }

    template<typename Iterator>
    AoSVector(Iterator start, Iterator end)
        : AoSVector{check_distance_size_t(start, end)}
    {
        std::copy(start, end, particles_.begin());
    }

    AoSVector(AoSVector const& from) = default;
    AoSVector(AoSVector&& from)      = default;
    AoSVector& operator=(AoSVector&& from) = default;
    AoSVector& operator=(AoSVector const& from) = default;

    template<typename T>
    struct iterator_impl;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    auto size() const { return particles_.size(); }

    void clear() { particles_.clear(); }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }

    auto& operator[](std::size_t i) const { return particles_.data()[i]; }
    auto& operator[](std::size_t i) { return particles_.data()[i]; }

    bool operator==(This const& that) const { return (this->particles_ == that.particles_); }

    auto begin() const { return const_iterator{particles_.begin(), *this}; }
    auto begin() { return iterator{particles_.begin(), *this}; }

    auto end() const { return const_iterator{particles_.end(), *this}; }
    auto end() { return iterator{particles_.end(), *this}; }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

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
        return particles_.erase(particles_.begin() + first.idx(), particles_.begin() + last.idx());
    }

    Particle_t& emplace_back() { return particles_.emplace_back(); }



    Particle_t& emplace_back(Particle_t&& p)
    {
        return particles_.emplace_back(std::forward<Particle_t>(p));
    }

    Particle_t& emplace_back(Particle_t const& p)
    {
        return particles_.emplace_back(std::forward<Particle_t>(p));
    }

    template<typename... Args>
    Particle_t& emplace_back(Args const&... args)
    {
        return particles_.emplace_back(args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(Args&&... args)
    {
        return particles_.emplace_back(Particle_t{args...});
    }

    void push_back(Particle_t const& p) { particles_.push_back(p); }

    void push_back(Particle_t&& p) { particles_.push_back(std::forward<Particle_t>(p)); }

    // template<typename That>
    void swap(This& that) { std::swap(this->particles_, that.particles_); }




protected:
    container_type particles_;
};



template<typename Super_>
struct AoSParticles : public Super_
{
    using Super = Super_;
    using This  = AoSParticles<Super>;

    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = Super::dimension;

    template<typename S = Super_>
    static constexpr bool is_vector()
    {
        return std::is_same_v<S, AoSVector<dimension>>;
    }

    using Super::particles_;
    using container_type = typename Super::container_type;
    using Particle_t     = typename Super::Particle_t;

    template<std::size_t size>
    using array_type  = AoSParticles<AoSArray<dimension, size>>;
    using vector_type = AoSParticles<AoSVector<dimension>>;



    AoSParticles() {}

    template<typename S = Super, std::enable_if_t<is_vector<S>(), bool> = true>
    AoSParticles(std::size_t size)
        : Super(size)
    {
    }

    template<typename Particle_t, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(std::size_t size, Particle_t&& particle)
        : Super(size, std::forward<Particle_t>(particle))
    {
    }


    template<typename Iterator, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(Iterator start, Iterator end)
        : Super{start, end}
    {
    }

    AoSParticles(AoSParticles const& from) = default;
    AoSParticles(AoSParticles&& from)      = default;
    AoSParticles& operator=(AoSParticles&& from) = default;
    AoSParticles& operator=(AoSParticles const& from) = default;


    template<typename T, bool is_const = false> // :(
    auto static constexpr iterator_type()
    {
        throw std::runtime_error("never to be called in non constexpr context");
        if constexpr (is_vector<>())
            return std::conditional_t<is_const, typename Super::template iterator_impl<T const>*,
                                      typename Super::template iterator_impl<T>*>{nullptr};
        else
            return std::conditional_t<is_const, typename Super::const_iterator,
                                      typename Super::iterator>{nullptr};
    }

    using iterator       = std::decay_t<decltype(*iterator_type<This>())>;
    using const_iterator = std::decay_t<decltype(*iterator_type<This, true>())>;

    template<typename It, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    auto erase(It a, It b)
    {
        return Super::erase(a, b);
    }

    auto begin() const
    {
        if constexpr (is_vector<>())
            return const_iterator{particles_.begin(), *this};
        else
            return Super::begin();
    }
    auto begin()
    {
        if constexpr (is_vector<>())
            return iterator{particles_.begin(), *this};
        else
            return Super::begin();
    }

    auto end() const
    {
        if constexpr (is_vector<>())
            return const_iterator{particles_.end(), *this};
        else
            return Super::end();
    }
    auto end()
    {
        if constexpr (is_vector<>())
            return iterator{particles_.end(), *this};
        else
            return Super::end();
    }

    auto& weight(std::size_t i) const { return particles_[i].weight(); }
    auto& weight(std::size_t i) { return particles_[i].weight(); }

    auto& charge(std::size_t i) const { return particles_[i].charge(); }
    auto& charge(std::size_t i) { return particles_[i].charge(); }

    auto& iCell(std::size_t i) const { return particles_[i].iCell(); }
    auto& iCell(std::size_t i) { return particles_[i].iCell(); }

    auto& delta(std::size_t i) const { return particles_[i].delta(); }
    auto& delta(std::size_t i) { return particles_[i].delta(); }

    auto& v(std::size_t i) const { return particles_[i].v(); }
    auto& v(std::size_t i) { return particles_[i].v(); }

    auto& E(std::size_t i) const { return particles_[i].E_; }
    auto& E(std::size_t i) { return particles_[i].E_; }

    auto& B(std::size_t i) const { return particles_[i].B_; }
    auto& B(std::size_t i) { return particles_[i].B_; }

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }


    template<typename S = Super, std::enable_if_t<is_vector<S>(), bool> = 0>
    auto capacity() const
    {
        return particles_.capacity();
    }
    auto size() const { return particles_.size(); }

    template<typename IndexRange, typename Predicate>
    auto partition(IndexRange&& range, Predicate&& pred)
    {
        return std::partition(range.begin(), range.end(), pred);
    }

    auto& vector() { return particles_; }
    auto& vector() const { return particles_; }


    auto constexpr static size_of_particle() { return sizeof(Particle_t); }
};


template<std::size_t dim, std::size_t size>
using AoSArrayParticles = AoSParticles<AoSArray<dim, size>>;

template<std::size_t dim>
using AoSVectorParticles = AoSParticles<AoSVector<dim>>;

template<typename Super_>
class AoSMappedParticles : public Super_
{
    using Super = Super_;
    using This  = AoSMappedParticles<Super>;
    using Super::particles_;

    static_assert(Super::template is_vector<>(), "not supported otherwise");

public:
    static constexpr auto dimension = Super::dimension;

    using box_t      = Box<int, dimension>;
    using CellMap_t  = CellMap<dimension, int>;
    using Particle_t = typename Super::Particle_t;

    template<std::size_t size>
    using array_type  = typename Super::template array_type<size>;
    using vector_type = typename Super::vector_type;


    AoSMappedParticles(box_t box)
        : Super_{}
        , box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }

    AoSMappedParticles(box_t box, std::size_t size)
        : Super_(size)
        , box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }


    template<typename T, bool is_const = false> // :(
    auto static constexpr iterator_type()
    {
        throw std::runtime_error("never to be called in non constexpr context");
        return std::conditional_t<is_const, typename Super::template iterator_impl<T const>*,
                                  typename Super::template iterator_impl<T>*>{nullptr};
    }

    using iterator       = std::decay_t<decltype(*iterator_type<This>())>;
    using const_iterator = std::decay_t<decltype(*iterator_type<This, true>())>;

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


    auto nbr_particles_in(box_t const& box) const { return cellMap_.size(box); }

    void export_particles(box_t const& box, This& dest) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles");
        cellMap_.export_to(box, *this, dest);
    }

    template<typename Dest, typename Fn>
    void export_particles(box_t const& box, Dest& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn)");
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (box, vector, Fn)");
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(This& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn,vector)");
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

    auto begin() const { return const_iterator{particles_.begin(), *this}; }
    auto begin() { return iterator{particles_.begin(), *this}; }

    auto end() const { return const_iterator{particles_.end(), *this}; }
    auto end() { return iterator{particles_.end(), *this}; }

protected:
    box_t box_{};
    mutable CellMap_t cellMap_;
};

template<std::size_t dim>
using AoSMappedVectorParticles = AoSMappedParticles<AoSParticles<AoSVector<dim>>>;


template<std::size_t dim>
template<typename T>
struct AoSVector<dim>::iterator_impl
    : public wrapped_iterator<T, typename std::decay_t<T>::container_type>
{
    static constexpr auto is_contiguous = false;
    static constexpr auto dimension     = dim;

    using Super = wrapped_iterator<T, typename std::decay_t<T>::container_type>;

    iterator_impl() {}

    template<typename Iterator>
    iterator_impl(Iterator iter, T& container)
        : Super{{iter}, &container}
    {
    }
    iterator_impl(iterator_impl const& that) = default;

    iterator_impl& operator=(iterator_impl const& other) = default;

    bool operator==(iterator_impl const& other) const
    {
        bool superbool = (static_cast<Super const&>(*this) == static_cast<Super const&>(other));
        return superbool;
    }

    auto operator+(std::size_t i)
    {
        iterator_impl copy = *this;
        static_cast<Super&>(copy) += i;
        return copy;
    }

    auto& weight() { return (*this)->weight(); }
    auto& weight() const { return (*this)->weight(); }

    auto& charge() { return (*this)->charge(); }
    auto& charge() const { return (*this)->charge(); }

    auto& iCell() { return (*this)->iCell(); }
    auto& iCell() const { return (*this)->iCell(); }

    auto& delta() { return (*this)->delta(); }
    auto& delta() const { return (*this)->delta(); }

    auto& v() { return (*this)->v(); }
    auto& v() const { return (*this)->v(); }

    template<typename Weight>
    void weight(Weight const& weight)
    {
        (*this)->weight_ = weight;
    }

    template<typename Charge>
    void charge(Charge const& charge)
    {
        (*this)->charge_ = charge;
    }

    template<typename ICell>
    void iCell(ICell const& iCell)
    {
        (*this)->iCell_ = iCell;
    }

    template<typename Delta>
    void delta(Delta const& delta)
    {
        (*this)->delta_ = delta;
    }

    template<typename V>
    void v(V const& v)
    {
        (*this)->v_ = v;
    }
};

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP */
