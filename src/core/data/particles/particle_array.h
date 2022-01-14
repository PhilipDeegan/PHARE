#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <utility>
#include <vector>

#include "particle.h"
#include "core/utilities/point/point.h"
#include "core/utilities/cellmap.h"
#include "core/utilities/iterators.h"
#include "core/logger.h"

namespace PHARE::core
{
template<std::size_t dim, std::size_t size>
struct AoSArray
{
    using container_type = std::array<Particle<dim>, size>;

    container_type particles_;
};

template<std::size_t dim>
struct AoSVector
{
    using container_type = std::vector<Particle<dim>>;
    using value_type     = typename container_type::value_type;

    AoSVector() {}

    AoSVector(std::size_t size)
        : particles_(size)
    {
    }

    AoSVector(std::size_t size, value_type&& particle)
        : particles_(size, particle)
    {
    }

    container_type particles_;
};


template<std::size_t dim, typename Super_ = AoSVector<dim>>
struct AoSParticles : public Super_
{
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;

    using This  = AoSParticles<dim, Super_>;
    using Super = Super_;
    using Super::particles_;

    template<typename T>
    auto static constexpr is_vector()
    {
        return std::is_same_v<T, AoSVector<dim>>;
    }

    AoSParticles() {}

    template<typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(std::size_t size)
        : Super(size)
    {
    }

    template<typename Particle_t, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(std::size_t size, Particle_t&& particle)
        : Super(size, std::forward<Particle_t>(particle))
    {
    }

    template<typename T>
    struct iterator_t;
    using iterator       = iterator_t<This>;
    using const_iterator = iterator_t<This const>;
};

template<std::size_t dim, typename OuterSuper>
template<typename T>
struct AoSParticles<dim, OuterSuper>::iterator_t
    : wrapped_iterator<T, typename OuterSuper::container_type>
{
    static constexpr auto is_contiguous = false;
    static constexpr auto dimension     = dim;

    using Super = wrapped_iterator<T, typename OuterSuper::container_type>;

    template<typename Iterator>
    iterator_t(Iterator iter, T* container)
        : Super{iter, container}
    {
    }
    iterator_t(iterator_t const& that) = default;

    iterator_t& operator=(iterator_t const& other) = default;

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
    void weight(Weight weight)
    {
        (*this)->weight_ = weight;
    }

    template<typename Charge>
    void charge(Charge charge)
    {
        (*this)->charge_ = charge;
    }

    template<typename ICell>
    void iCell(ICell iCell)
    {
        (*this)->iCell_ = iCell;
    }

    template<typename Delta>
    void delta(Delta delta)
    {
        (*this)->delta_ = delta;
    }

    template<typename V>
    void v(V v)
    {
        (*this)->v_ = v;
    }
};


template<std::size_t dim>
class ParticleArray : public AoSParticles<dim, AoSVector<dim>>
{
public:
    static constexpr bool is_contiguous              = false;
    static constexpr auto dimension                  = dim;
    static constexpr std::size_t cellmap_bucket_size = 100;
    using This                                       = ParticleArray<dim>;
    using Super                                      = AoSParticles<dim, AoSVector<dim>>;
    using Particle_t                                 = Particle<dim>;
    using Vector                                     = typename Super::container_type;
    using Super::particles_;

    template<std::size_t size>
    using array_type = AoSParticles<dim, AoSArray<dim, size>>;


private:
    using cell_map_t = CellMap<dim, Particle_t, cellmap_bucket_size, int, Point<int, dim>>;


public:
    using value_type = Particle_t;
    using box_t      = Box<int, dim>;

    template<typename T, typename cell_map_t>
    struct ParticleArrayIterator;

    using iterator       = ParticleArrayIterator<This, cell_map_t>;
    using const_iterator = ParticleArrayIterator<This const, cell_map_t>;


    ParticleArray() {}
    ParticleArray(std::size_t size)
        : Super{size}
    {
    }

    ParticleArray(std::size_t size, Particle_t&& particle)
        : Super{size, std::forward<Particle_t>(particle)}
    {
    }

    auto& iCell(std::size_t i) const { return particles_[i].iCell(); }
    auto& iCell(std::size_t i) { return particles_[i].iCell(); }


    std::size_t size() const { return particles_.size(); }
    std::size_t capacity() const { return particles_.capacity(); }

    void clear()
    {
        clean_ = false;
        return particles_.clear();
    }
    void reserve(std::size_t newSize)
    {
        clean_ = false;
        return particles_.reserve(newSize);
    }
    void resize(std::size_t newSize)
    {
        clean_ = false;
        return particles_.resize(newSize);
    }

    auto& operator[](std::size_t i) const { return particles_[i]; }
    auto& operator[](std::size_t i) { return particles_[i]; }

    bool operator==(ParticleArray<dim> const& that) const
    {
        return (this->particles_ == that.particles_);
    }

    auto begin() const { return const_iterator{particles_.begin(), *this, cell_map_}; }
    auto begin() { return iterator{particles_.begin(), *this, cell_map_}; }

    auto end() const { return const_iterator{particles_.end(), *this, cell_map_}; }
    auto end() { return iterator{particles_.end(), *this, cell_map_}; }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

    auto back()
    {
        clean_ = false;
        return particles_.back();
    }
    auto front()
    {
        clean_ = false;
        return particles_.front();
    }

    iterator erase(iterator position)
    {
        clean_ = false;
        return particles_.erase(position);
    }
    iterator erase(iterator first, iterator last)
    {
        clean_ = false;
        // return particles_.erase(first, last);
        return iterator{particles_.erase(first, last), this};
    }

    Particle_t& emplace_back()
    {
        clean_ = false;
        return particles_.emplace_back();
    }
    Particle_t& emplace_back(Particle_t&& p)
    {
        clean_ = false;
        return particles_.emplace_back(p);
    }
    Particle_t& emplace_back(Particle_t const& p)
    {
        clean_ = false;
        return particles_.emplace_back(p);
    }

    template<typename... Args>
    Particle_t& emplace_back(Args const&... args)
    {
        clean_ = false;
        return particles_.emplace_back(args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(Args&&... args)
    {
        clean_ = false;
        return particles_.emplace_back(Particle_t{args...});
    }

    void push_back(Particle_t const& p)
    {
        clean_ = false;
        particles_.push_back(p);
    }
    void push_back(Particle_t&& p)
    {
        clean_ = false;
        particles_.push_back(std::forward<Particle_t>(p));
    }

    void swap(ParticleArray<dim>& that)
    {
        clean_ = false;
        std::swap(this->particles_, that.particles_);
    }

    void map_particles() const
    {
        cell_map_.add(particles_);
        clean_ = true;
    }
    void empty_map() { cell_map_.empty(); }


    auto nbr_particles_in(box_t const& box) const
    {
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        auto s = cell_map_.size(box);
        return s;
    }


    auto select(box_t const& box) const
    {
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        return cell_map_.select(box);
    }


    void export_particles(box_t const& box, ParticleArray<dim>& dest) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest.particles_);
    }

    template<typename Fn>
    void export_particles(box_t const& box, ParticleArray<dim>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn)");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest.particles_, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn,vector)");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest, std::forward<Fn>(fn));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cell_map_.print(cell);
    }

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

    auto constexpr static size_of() { return sizeof(Particle_t); }

private:
    bool mutable clean_{false};
    mutable cell_map_t cell_map_;
};


template<std::size_t dim>
template<typename T, typename cell_map_t>
struct ParticleArray<dim>::ParticleArrayIterator
    : public std::decay_t<T>::Super::template iterator_t<T>
{
    static constexpr auto is_contiguous = false;
    static constexpr auto dimension     = dim;

    using Super = typename std::decay_t<T>::Super::template iterator_t<T>;

    template<typename Iterator>
    ParticleArrayIterator(Iterator iter, T& container, cell_map_t& cellmap)
        : Super{iter, &container}
        , cm{&cellmap}
    {
    }
    ParticleArrayIterator(ParticleArrayIterator const& that) = default;

    template<typename Cell>
    void change_icell(Cell const& newCell)
    {
        auto particleIndex = std::distance(this->container->begin(), *this);
        cm->update(*(this->container), particleIndex, newCell);
        (*this)->iCell = newCell;
    }

    bool operator==(ParticleArrayIterator const& other) const
    {
        bool superbool = (static_cast<Super const&>(*this) == static_cast<Super const&>(other));
        bool cmb       = (cm == other.cm);
        return superbool and cmb;
    }
    ParticleArrayIterator& operator=(ParticleArrayIterator const& other) = default;

    cell_map_t* cm;
};

} // namespace PHARE::core


namespace PHARE
{
namespace core
{
    template<std::size_t dim>
    void empty(ParticleArray<dim>& array)
    {
        array.clear();
    }

    template<std::size_t dim>
    void swap(ParticleArray<dim>& array1, ParticleArray<dim>& array2)
    {
        array1.swap(array2);
    }

} // namespace core
} // namespace PHARE



#endif
