#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include "core/def.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"
#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/particle_array_selector.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"

#include "core/utilities/equality.hpp"

namespace PHARE::core
{
/*
template<std::size_t dim,                         //
         auto layout_mode  = LayoutMode::AoS,     //
         auto storage_mode = StorageMode::VECTOR, //
         auto alloc_mode   = AllocatorMode::CPU   // only matters for vector
         >
class ParticleArrayInternals
*/

template<std::size_t dim, typename internals>
class ParticleArray : public internals::value_type
{
public:
    using Super                        = typename internals::value_type;
    using value_type                   = typename ParticleDefaults<dim>::Particle_t;
    auto static constexpr dimension    = dim;
    auto static constexpr alloc_mode   = internals::alloc_mode;
    auto static constexpr is_mapped    = internals::is_mapped; // torm
    auto static constexpr layout_mode  = internals::layout_mode;
    auto static constexpr storage_mode = internals::storage_mode;
    auto static constexpr type_id      = internals::type_id;
    auto static constexpr impl         = internals::impl;
    using view_t                       = ParticleArray<
                              dim, ParticleArrayInternals<dim, layout_mode, StorageMode::SPAN, alloc_mode, impl>>;

    std::string static id()
    {
        std::stringstream ss;
        ss << magic_enum::enum_name(layout_mode) << "," << magic_enum::enum_name(alloc_mode) << ","
           << magic_enum::enum_name(storage_mode);
        return ss.str();
    }

    template<typename... Args>
    ParticleArray(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        return Super::partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    auto& replace_from(ParticleArray const& that)
    {
        if (this != &that)
            Super::replace_from(that);
        return *this;
    }

    auto view() { return view_t{*this}; }
    auto view() const { return view_t{*this}; }

    auto view(std::size_t i) // to take only i particles and ignore the rest
    {
        view_t v{*this};
        v.super().resize(i);
        return v;
    }

    auto view(std::size_t const start, std::size_t const size)
    {
        return view_t{*this, start, size};
    }

    auto operator*() { return view(); }
    auto operator*() const { return view(); }

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super& super() const { return *this; }

    template<std::size_t _dim, typename _internals>
    friend std::ostream& operator<<(std::ostream& out, ParticleArray<dim, _internals> const&);
};

template<std::size_t dim>
using AoSParticleArray = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoS>>;

template<std::size_t dim>
using AoSMappedParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoSMapped>>;

template<std::size_t dim>
using SoAParticleArray = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::SoA>>;


template<std::size_t dim, typename internals>
std::ostream& operator<<(std::ostream& out, ParticleArray<dim, internals> const& arr)
{
    if constexpr (ParticleArray<dim, internals>::layout_mode == LayoutMode::SoAPC) {}
    else
        for (auto const& p : arr)
            out << p.copy();
    return out;
}

template<typename Particles, typename T, std::size_t S, typename Particle_t = Particle<S>>
auto make_particles(Box<T, S> const& box, std::size_t const s = 0, Particle_t const& particle = {})
{
    using enum LayoutMode;
    if constexpr (Particles::is_mapped)
        return Particles{box, s, particle};
    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC))
        return Particles{box};
    else
        return make_particles<Particles>(s, particle);
}

template<typename Particles, typename Particle_t = Particle<Particles::dimension>>
auto make_particles(std::size_t const s = 0, Particle_t const& particle = {})
{
    return Particles{s, particle};
}


template<typename Particles, typename GridLayout_t>
auto make_particles(GridLayout_t const& layout)
{
    auto constexpr static ghost_cells = GridLayout_t::nbrParticleGhosts();

    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC, AoSTS, SoATS, SoAVXTS))
        return Particles{layout.AMRBox(), ghost_cells};

    else if constexpr (Particles::is_mapped)
        return Particles{grow(layout.AMRBox(), ghost_cells + 1)};

    else
        return Particles{};
}

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

template<bool binary_equal = false, typename P0, typename P1>
EqualityReport particle_compare(P0 const& p0, P1 const& p1, std::size_t const i = 0)
{
    std::string idx = std::to_string(i);
    if (p0.iCell() != p1.iCell())
        return EqualityReport{false, "icell mismatch at index: " + idx, i};
    if constexpr (binary_equal)
    {
        if (p0.delta() != p1.delta())
            return EqualityReport{false, "delta mismatch at index: " + idx, i};
        if (p0.v() != p1.v())
            return EqualityReport{false, "v mismatch at index: " + idx, i};
    }
    else
    {
        if (!float_equals(p0.delta(), p1.delta()))
            return EqualityReport{false, "delta mismatch at index: " + idx, i};
        if (!float_equals(p0.v(), p1.v()))
            return EqualityReport{false, "v mismatch at index: " + idx, i};
    }

    return EqualityReport{true};
}


template<bool binary_equal = false, typename P0, typename P1>
EqualityReport compare_particles(P0 const& ref, P1 const& cmp)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ref.size()) + " vs "
                                         + std::to_string(cmp.size())};

    auto rit      = ref.begin();
    auto cit      = cmp.begin();
    std::size_t i = 0;

    for (; rit != ref.end(); ++rit, ++cit, ++i)
        if (auto const eq = particle_compare<binary_equal>(*rit, *cit, i); !eq)
            return eq;

    return EqualityReport{true};
}

template<std::size_t D, typename I0, typename I1>
EqualityReport operator==(ParticleArray<D, I0> const& p0, ParticleArray<D, I1> const& p1)
{
    auto report = compare_particles(p0, p1);
    if (!report)
    {
        PHARE_LOG_LINE_STR(p0[report.idx].copy());
        PHARE_LOG_LINE_STR(p1[report.idx].copy());
    }
    return report;
}

template<std::size_t D, typename I0>
EqualityReport operator==(ParticleArray<D, I0> const& p0, std::vector<Particle<D>> const& p1)
{
    return compare_particles(p0, p1);
}

struct Sortings
{ // to be used in contexpr fashion

    bool by_delta = true;
};

template<Sortings S = Sortings{}, std::size_t D, typename I0, typename Box_t>
auto& sort(ParticleArray<D, I0>& ps, Box_t const& box)
{
    std::string_view constexpr static FN_ID = "sort,";
    auto constexpr function_id = join_string_views_v<FN_ID, ParticleArray<D, I0>::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    if (ps.size() == 0)
        return ps;

    ParticleArraySorter<ParticleArray<D, I0>> sorter{ps, box};
    sorter();
    if constexpr (S.by_delta)
        sorter.by_delta();

    using enum LayoutMode;
    if constexpr (any_in(ParticleArray<D, I0>::layout_mode, AoSPC, SoAPC))
        ps.reset_index_wrapper_map();

    return ps;
}


} // namespace PHARE::core

#endif
