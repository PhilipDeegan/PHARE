#ifndef PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_CPU_HPP
#define PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_CPU_HPP

// no includes, is included

#include <stdexcept>

namespace PHARE::core::detail
{

template<typename V>
bool constexpr el_wise_less(V const& v0, V const& v1)
{
    for (std::int16_t i = v0.size() - 1; i >= 0; --i)
    {
        if (v0[i] < v1[i])
            return true;
        if (v0[i] > v1[i])
            return false;
    }
    return false;
}

template<typename V>
bool constexpr el_wise_gr8r(V const& v0, V const& v1)
{
    for (std::int16_t i = v0.size() - 1; i >= 0; --i)
    {
        if (v0[i] > v1[i])
            return true;
        if (v0[i] < v1[i])
            return false;
    }
    return false;
}

template<typename ParticleArray>
class ParticleSorter<AllocatorMode::CPU, /*impl = */ 0, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    ParticleSorter(ParticleArray& particles_, box_t domain)
        : particles{particles_}
        , domain_box{domain}
    {
    }



    ParticleSorter& operator()(std::int64_t const& l, std::int64_t const& r) // basically quicksort
    {
        using enum LayoutMode;

        if constexpr (ParticleArray::layout_mode == SoA)
        {
            PHARE_WITH_THRUST_ELSE_THROW(                      //
                auto it = SoAIteratorAdaptor::make(particles); //
                std::sort(it + l, it + r,
                          [cf = cell_flattener](auto const& a, auto const& b) {
                              return cf(SoAIteratorAdaptor::iCell(a))
                                     < cf(SoAIteratorAdaptor::iCell(b));
                          }); //
            )
        }
        else if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS, SoAVXTS))
        {
            // fix soa - not finished
            for (auto& tile : particles())
            {
                std::sort(tile().begin(), tile().end(),
                          [cf = cell_flattener](auto const& a, auto const& b) {
                              return cf(a.iCell()) < cf(b.iCell());
                          });
            }
        }
        else
        {
            std::sort(particles.begin() + l, particles.begin() + r,
                      [cf = cell_flattener](auto const& a, auto const& b) {
                          return cf(a.iCell()) < cf(b.iCell());
                      });
        }

        return *this;
    }


    ParticleSorter& operator()()
    {
        (*this)(0, particles.size());
        return *this;
    }

    auto static constexpr by_deltas()
    {
        return [](auto const& a, auto const& b) -> bool {
            return as_tuple(a.delta()) < as_tuple(b.delta());
        };
    }

    void by_deltas(std::uint64_t const& l, std::uint64_t const& r)
    {
        using enum LayoutMode;

        if constexpr (ParticleArray::layout_mode == SoA)
        {
            PHARE_WITH_THRUST_ELSE_THROW( //
                auto it = SoAIteratorAdaptor::make(particles);
                std::sort(it + l, it + r,
                          [](auto const& a, auto const& b) -> bool {
                              return as_tuple(SoAIteratorAdaptor::delta(a))
                                     < as_tuple(SoAIteratorAdaptor::delta(b));
                          }); //
            )
        }
        else if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS, SoAVXTS))
        {
            throw std::runtime_error("no sort for AoSTS/SoATS");
        }
        else
        {
            std::sort(particles.begin() + l, particles.begin() + r, by_deltas());
        }
    }


    ParticleSorter& by_delta()
    {
        // assumes already sorted by icell
        if (particles.size() == 0)
            return *this;

        using enum LayoutMode;

        if constexpr (any_in(ParticleArray::layout_mode, AoSPC, SoAPC))
            for (auto const& bix : particles.local_box())
                std::sort(particles(bix).begin(), particles(bix).end(), by_deltas());
        else if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS, SoAVXTS))
        {
            // fix soa - not finished
            for (auto& tile : particles())
                std::sort(tile().begin(), tile().end(), by_deltas());
        }
        else
        {
            auto const end = particles.end();
            auto beg       = particles.begin();
            auto lst       = particles.begin();
            while (lst != end)
            {
                lst = beg + 1;
                while (lst != end and lst.iCell() == beg.iCell())
                    ++lst;
                by_deltas(it_dist(particles.begin(), beg), it_dist(particles.begin(), lst));
                beg = lst;
            }
        }

        return *this;
    }


private:
    ParticleArray& particles;
    box_t domain_box;
    CellFlattener<box_t> cell_flattener{domain_box};


    template<typename V>
    bool constexpr flat_less(V const& v0, V const& v1) const
    {
        return cell_flattener(v0) < cell_flattener(v1);
    }

    template<typename V>
    bool constexpr flat_gr8r(V const& v0, V const& v1) const
    {
        return cell_flattener(v0) > cell_flattener(v1);
    }
};


template<typename ParticleArray>
class ParticleSorter<AllocatorMode::CPU, /*impl = */ 1, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    void operator()(std::int64_t const& l, std::int64_t const& r) // basically quicksort
    {
        auto i          = l;
        auto j          = r;
        auto const half = particles.iCell((l + r) / 2);
        do
        {
            while (el_wise_less(particles.iCell(i), half))
                i++;
            while (el_wise_gr8r(particles.iCell(j), half))
                j--;
            if (i <= j)
            {
                particles.swap(i, j);
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


    ParticleArray& particles;
    box_t domain_box;
};



} // namespace PHARE::core::detail

#endif /*PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_CPU_HPP*/
