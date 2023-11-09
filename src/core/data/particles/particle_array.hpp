#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP


#include <cstddef>
#include <utility>
#include <vector>

#include "core/utilities/indexer.hpp"
#include "particle.hpp"
#include "particle_array_defs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/def.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;
    using This                          = ParticleArray<dim>;
    using Particle_t                    = Particle<dim>;
    using Vector                        = std::vector<Particle_t>;

private:
    using CellMap_t   = CellMap<dim, int>;
    using IndexRange_ = IndexRange<This>;


public:
    using value_type     = Particle_t;
    using box_t          = Box<int, dim>;
    using iterator       = typename Vector::iterator;
    using const_iterator = typename Vector::const_iterator;



public:
    ParticleArray(box_t box)
        : box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }

    ParticleArray(box_t box, std::size_t size)
        : particles_(size)
        , box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }

    ParticleArray(ParticleArray const& from)            = default;
    ParticleArray(ParticleArray&& from)                 = default;
    ParticleArray& operator=(ParticleArray&& from)      = default;
    ParticleArray& operator=(ParticleArray const& from) = default;

    NO_DISCARD std::size_t size() const { return particles_.size(); }
    NO_DISCARD std::size_t capacity() const { return particles_.capacity(); }

    void clear()
    {
        particles_.clear();
        cellMap_.clear();
    }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }

    NO_DISCARD auto const& operator[](std::size_t i) const { return particles_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) { return particles_[i]; }

    NO_DISCARD bool operator==(ParticleArray<dim> const& that) const
    {
        return (this->particles_ == that.particles_);
    }

    NO_DISCARD auto begin() const { return particles_.begin(); }
    NO_DISCARD auto begin() { return particles_.begin(); }

    NO_DISCARD auto end() const { return particles_.end(); }
    NO_DISCARD auto end() { return particles_.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

    NO_DISCARD auto back() { return particles_.back(); }
    NO_DISCARD auto front() { return particles_.front(); }

    auto erase(IndexRange_& range) { cellMap_.erase(particles_, range); }
    auto erase(IndexRange_&& range)
    {
        // TODO move ctor for range?
        cellMap_.erase(std::forward<IndexRange_>(range));
    }

    iterator erase(iterator first, iterator last)
    {
        // should we erase particles indexes associated with these iterators from the cellmap?
        // probably it does not matter if not. The reason is that
        // particles erased from the particlearray are so because they left
        // the patch cells to an outside cell.
        // But in principle that cell will never be accessed because it is outside the patch.
        // The only thing "bad" if these indexes are not deleted is that the
        // size of the cellmap becomes unequal to the size of the particleArray.
        // but  ¯\_(ツ)_/¯
        return particles_.erase(first, last);
    }
    iterator erase(const_iterator first, const_iterator last)
    {
        return particles_.erase(first, last);
    }


    Particle_t& emplace_back()
    {
        auto& part = particles_.emplace_back();
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }



    Particle_t& emplace_back(Particle_t&& p)
    {
        auto& part = particles_.emplace_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    void push_back(Particle_t const& p)
    {
        particles_.push_back(p);
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void push_back(Particle_t&& p)
    {
        particles_.push_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void swap(ParticleArray<dim>& that) { std::swap(this->particles_, that.particles_); }

    void map_particles() const { cellMap_.add(particles_); }
    void empty_map() { cellMap_.empty(); }


    NO_DISCARD auto nbr_particles_in(box_t const& box) const { return cellMap_.size(box); }

    void export_particles(box_t const& box, ParticleArray<dim>& dest) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles");
        cellMap_.export_to(box, particles_, dest);
    }

    template<typename Fn>
    void export_particles(box_t const& box, ParticleArray<dim>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (box, vector, Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(This& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn,vector)");
        cellMap_.export_if(particles_.data(), dest, std::forward<Predicate>(pred));
    }


    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t particleIndex)
    {
        auto oldCell                    = particles_[particleIndex].iCell;
        particles_[particleIndex].iCell = newCell;
        if (!box_.isEmpty())
        {
            cellMap_.update(particles_, particleIndex, oldCell);
        }
    }


    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        return cellMap_.partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cellMap_.print(cell);
    }


    NO_DISCARD bool is_mapped() const
    {
        bool ok = true;
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

    void sortMapping() const { cellMap_.sort(); }

    NO_DISCARD auto& vector() { return particles_; }
    NO_DISCARD auto& vector() const { return particles_; }

    auto& box() { return box_; }
    auto& box() const { return box_; }

    auto& ppc() const { return ppc_; }
    auto& ppc() { return ppc_; }

    template<typename ICell>
    auto& ppc(ICell const& icell)
    {
        return ppc_(toLocal(icell));
    }
    template<typename ICell>
    auto& ppc(ICell const& icell) const
    {
        return ppc_(toLocal(icell));
    }

    template<typename ICell>
    auto& ppc_offsets(ICell const& icell)
    {
        return ppc_offsets_(toLocal(icell));
    }
    template<typename ICell>
    auto& ppc_offsets(ICell const& icell) const
    {
        return ppc_offsets_(toLocal(icell));
    }

    auto& update_from(This const& that)
    {
        this->resize(that.size());
        std::copy(that.begin(), that.end(), this->begin());
        this->box_ = that.box_;
        // this->cellMap_ = that.cellMap_;
        reset_sorting_bits();
        return *this;
    }
    auto& set_as(This const& that)
    {
        // this->resize(that.size());
        // std::copy(that.begin(), that.end(), this->begin());
        this->box_ = that.box_;
        // this->cellMap_ = that.cellMap_;
        reset_sorting_bits();
        return *this;
    }

private:
    void reset_sorting_bits()
    {
        auto shape = this->box_.shape().template toArray<std::uint32_t>();
        ppc_.update_from(shape);
        ppc_offsets_.update_from(shape);
    }

    // updater can put domain particles in first ghost cell for domain particles
    std::int32_t constexpr static ghostish_size = 0;

    template<typename ICell>
    auto toLocal(ICell icell) const
    {
        for_N<dim>([&](auto ic) {
            constexpr auto i = ic();
            icell[i] -= (box_.lower[i] - ghostish_size);
            DEBUG_ABORT_IF(icell[i] < 0);
        });
        return icell;
    }

    Vector particles_;
    box_t box_;
    mutable CellMap_t cellMap_;

    NdArrayVector<dim, std::uint32_t> ppc_{box().shape().template toArray<std::uint32_t>()};
    NdArrayVector<dim, std::uint32_t> ppc_offsets_{box().shape().template toArray<std::uint32_t>()};
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


    template<std::size_t dim, bool OwnedState = true>
    struct ContiguousParticles
    {
        static constexpr bool is_contiguous    = true;
        static constexpr std::size_t dimension = dim;
        using ContiguousParticles_             = ContiguousParticles<dim, OwnedState>;

        template<typename T>
        using container_t = std::conditional_t<OwnedState, std::vector<T>, Span<T>>;

        template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
        ContiguousParticles(std::size_t s)
            : iCell(s * dim)
            , delta(s * dim)
            , weight(s)
            , charge(s)
            , v(s * 3)
        {
        }

        template<typename Container_int, typename Container_double>
        ContiguousParticles(Container_int&& _iCell, Container_double&& _delta,
                            Container_double&& _weight, Container_double&& _charge,
                            Container_double&& _v)
            : iCell{_iCell}
            , delta{_delta}
            , weight{_weight}
            , charge{_charge}
            , v{_v}
        {
        }

        NO_DISCARD std::size_t size() const { return weight.size(); }

        template<std::size_t S, typename T>
        NO_DISCARD static std::array<T, S>* _array_cast(T const* array)
        {
            return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
        }

        template<typename Return>
        NO_DISCARD Return _to(std::size_t i)
        {
            return {
                *const_cast<double*>(weight.data() + i),     //
                *const_cast<double*>(charge.data() + i),     //
                *_array_cast<dim>(iCell.data() + (dim * i)), //
                *_array_cast<dim>(delta.data() + (dim * i)), //
                *_array_cast<3>(v.data() + (3 * i)),
            };
        }

        NO_DISCARD auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
        NO_DISCARD auto view(std::size_t i) { return _to<ParticleView<dim>>(i); }

        NO_DISCARD auto operator[](std::size_t i) const { return view(i); }
        NO_DISCARD auto operator[](std::size_t i) { return view(i); }

        struct iterator
        {
            iterator(ContiguousParticles_* particles)
            {
                for (std::size_t i = 0; i < particles->size(); i++)
                    views.emplace_back((*particles)[i]);
            }

            iterator& operator++()
            {
                ++curr_pos;
                return *this;
            }

            NO_DISCARD bool operator!=(iterator const& other) const
            {
                return curr_pos != views.size();
            }
            NO_DISCARD auto& operator*() { return views[curr_pos]; }
            NO_DISCARD auto& operator*() const { return views[curr_pos]; }

            std::size_t curr_pos = 0;
            std::vector<ParticleView<dim>> views;
        };

        NO_DISCARD auto as_tuple()
        {
            return std::forward_as_tuple(weight, charge, iCell, delta, v);
        }
        NO_DISCARD auto as_tuple() const
        {
            return std::forward_as_tuple(weight, charge, iCell, delta, v);
        }

        NO_DISCARD auto begin() { return iterator(this); }
        NO_DISCARD auto cbegin() const { return iterator(this); }

        NO_DISCARD auto end() { return iterator(this); }
        NO_DISCARD auto cend() const { return iterator(this); }

        container_t<int> iCell;
        container_t<double> delta;
        container_t<double> weight, charge, v;
    };


    template<std::size_t dim>
    using ContiguousParticlesView = ContiguousParticles<dim, /*OwnedState=*/false>;

} // namespace core
} // namespace PHARE


namespace PHARE::core
{

template<typename ParticleArray>
class ParticleSorter
{
public:
    using box_t = typename ParticleArray::box_t;

    ParticleSorter(ParticleArray& particles_)
        : particles{particles_}
    {
    }

    void operator()(std::int64_t const& l, std::int64_t const& r) // basically quicksort
    {
        auto i          = l;
        auto j          = r;
        auto const half = particles[(l + r) / 2].iCell;
        do
        {
            while (el_wise_less(particles[i].iCell, half))
                i++;
            while (el_wise_gr8r(particles[j].iCell, half))
                j--;
            if (i <= j)
            {
                std::swap(particles[i], particles[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (l < j)
            (*this)(l, j);
        if (i < r)
            (*this)(i, r);
    }

    void operator()() { (*this)(0, particles.size() - 1); }

private:
    template<typename V>
    bool constexpr el_wise_less(V const& v0, V const& v1) const
    {
        return cell_flattener(v0) < cell_flattener(v1);
    }

    template<typename V>
    bool constexpr el_wise_gr8r(V const& v0, V const& v1) const
    {
        return cell_flattener(v0) > cell_flattener(v1);
    }

    ParticleArray& particles;
    CellFlattener<box_t> cell_flattener{particles.box()};
};




} // namespace PHARE::core


namespace std
{

template<std::size_t dim>
void sort_custom(PHARE::core::ParticleArray<dim>& particles)
{
    PHARE::core::ParticleSorter{particles}();
}

template<std::size_t dim>
auto& sort(PHARE::core::ParticleArray<dim>& particles)
{
    using box_t = typename PHARE::core::ParticleArray<dim>::box_t;
    PHARE::core::LocalisedCellFlattener<box_t> cell_flattener{grow(particles.box(), 1)};
    std::sort(particles.vector().begin(), particles.vector().end(),
              [&](auto const& a, auto const& b) {
                  return cell_flattener(a.iCell) < cell_flattener(b.iCell);
              });
    return particles;
}

template<std::size_t dim>
auto sort(PHARE::core::ParticleArray<dim>&& particles)
{
    return sort(particles);
}


} // namespace std


#endif
