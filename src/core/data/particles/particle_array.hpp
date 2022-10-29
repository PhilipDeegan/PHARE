#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include <cstddef>

#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"

#include "core/data/particles/particle_array_soa.hpp"
#include "core/data/particles/particle_array_aos.hpp"

namespace PHARE::core
{
template<std::size_t dim, bool aos_ = 1, bool mapped_ = 0>
class ParticleArrayInternals
{
    template<typename T>
    struct S
    {
        using value_type = T;
    };

    auto static constexpr resolve_t()
    {
        if constexpr (mapped)
        {
            static_assert(aos, "No mapped SOA version");
            return S<AoSMappedVectorParticles<dim>>{};
        }
        else
        {
            if constexpr (aos)
                return S<AoSVectorParticles<dim>>{};
            else
                return S<SoAVectorParticles<dim>>{};
        }
    }

public:
    static constexpr bool aos    = aos_;
    static constexpr auto mapped = mapped_;

    using value_type = typename decltype(resolve_t())::value_type;
};


template<std::size_t dim, typename internals = ParticleArrayInternals<dim>>
class ParticleArray : public internals::value_type
{
public:
    auto static constexpr is_mapped = internals::mapped;

    using This  = ParticleArray<dim, internals>;
    using Super = typename internals::value_type;
    using Super::erase;

    static constexpr bool is_contiguous = Super::is_contiguous;
    static constexpr auto dimension     = Super::dimension;

    using box_t       = Box<int, dimension>;
    using IndexRange_ = IndexRange<This>;

    template<std::size_t size>
    using array_type = typename Super::template array_type<size>;

    using const_iterator = typename Super::const_iterator;
    using iterator       = typename Super::iterator;
    using value_type     = typename Super::value_type;
    // using Super::export_particles;
    // using Super::nbr_particles_in;
    // using Super::sortMapping;


    template<typename I = internals, std::enable_if_t<!I::mapped, bool> = true>
    ParticleArray()
    {
    }

    template<typename I = internals, std::enable_if_t<!I::mapped, bool> = true>
    ParticleArray(std::size_t size)
        : Super{size}
    {
    }

    template<typename Particle_t, typename I = internals, std::enable_if_t<!I::mapped, bool> = true>
    ParticleArray(std::size_t size, Particle_t&& particle)
        : Super{size, std::forward<Particle_t>(particle)}
    {
    }

    template<typename It>
    ParticleArray(It start, It end)
        : Super{start, end}
    {
    }


    template<typename I = internals, std::enable_if_t<I::mapped, bool> = true>
    ParticleArray(box_t box)
        : Super{box}
    {
        assert(box.size() > 0);
    }

    template<typename I = internals, std::enable_if_t<I::mapped, bool> = true>
    ParticleArray(box_t box, std::size_t size)
        : Super{box, size}
    {
        assert(box.size() > 0);
    }

    ParticleArray(This const& from) = default;
    ParticleArray(This&& from)      = default;
    This& operator=(This&& from) = default;
    This& operator=(This const& from) = default;

    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        return Super::partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    auto erase(IndexRange_& range) { return Super::erase(range); }
    auto erase(IndexRange_&& range) { return Super::erase(std::forward<IndexRange_>(range)); }

    template<typename Iterator>
    auto erase(Iterator first, Iterator last)
    {
        return iterator{Super::erase(first, last), *this};
    }
};

template<std::size_t dim>
using AoSParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, /*aos =*/1, /*mapped_ =*/0>>;

template<std::size_t dim>
using AoSMappedParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, /*aos =*/1, /*mapped_ =*/1>>;

template<std::size_t dim>
using SoAParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, /*aos =*/0, /*mapped_ =*/0>>;


template<std::size_t dim, typename internals>
void empty(ParticleArray<dim, internals>& array)
{
    array.clear();
}


template<std::size_t dim, typename internals>
void swap(ParticleArray<dim, internals>& array1, ParticleArray<dim, internals>& array2)
{
    array1.swap(array2);
}


} // namespace PHARE::core




#endif
