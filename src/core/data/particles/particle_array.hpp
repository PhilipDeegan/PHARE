#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include "core/def.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"
#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/particle_array_selector.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"

#include "core/utilities/equality.hpp"
#include "core/utilities/monitoring.hpp"

#include <tuple>

namespace PHARE::core
{


template<auto opts, typename internals /* defaulted in details header */>
class ParticleArray : public internals::value_type
{
    using This = ParticleArray<opts, internals>;


    template<typename... Args>
    bool consteval static require()
    {
        using Tup = std::tuple<Args...>;

        bool constexpr base = std::is_constructible_v<Super, Args&&...>;
        if constexpr (std::tuple_size_v<Tup> > 0)
        {
            bool constexpr isself
                = std::is_same_v<std::decay_t<std::tuple_element_t<0, Tup>>, This>;
            // static_assert(!isself);
            return !isself and base;
        }
        else
            return base;
    }

public:
    using Super      = internals::value_type;
    using value_type = ParticleDefaults<opts.dim>::Particle_t;
    using view_t     = ParticleArray<opts.with_storage(StorageMode::SPAN)>;

    auto static constexpr options      = opts;
    auto static constexpr dimension    = opts.dim;
    auto static constexpr alloc_mode   = opts.alloc_mode;
    auto static constexpr layout_mode  = opts.layout_mode;
    auto static constexpr storage_mode = opts.storage_mode;
    auto static constexpr type_id      = internals::type_id;
    auto static constexpr is_mapped    = internals::is_mapped; // torm
    auto static inline mon             = MemoryMonitor{std::string{type_id}};

    std::string static id()
    {
        std::stringstream ss;
        ss << magic_enum::enum_name(layout_mode) << "," << magic_enum::enum_name(alloc_mode) << ","
           << magic_enum::enum_name(storage_mode);
        return ss.str();
    }

    ParticleArray(ParticleArray&& that)
        : Super{std::forward<Super>(that)}
    {
        mon.move();
    }
    ParticleArray(ParticleArray const& that)
        : Super{that}
    {
        mon.copy();
    }

    ParticleArray& operator=(ParticleArray&& that)
    {
        mon.move_assign();
        super() = std::move(that.super());
        return *this;
    }
    ParticleArray& operator=(ParticleArray const& that)
    {
        mon.copy_assign();
        super() = that.super();
        return *this;
    }

    template<typename... Args>
    ParticleArray(Args&&... args)
        requires(require<Args...>())
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
        mon.create();
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
    Super const& super() const { return *this; }

    template<auto _opts, typename _internals>
    friend std::ostream& operator<<(std::ostream& out, ParticleArray<_opts, _internals> const&);
};

template<std::size_t dim>
using AoSParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>;


template<std::size_t dim>
using AoSMappedParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoSMapped}>;

template<std::size_t dim>
using SoAParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::SoA}>;


template<auto opts, typename internals>
std::ostream& operator<<(std::ostream& out, ParticleArray<opts, internals> const& arr)
{
    if constexpr (ParticleArray<opts, internals>::layout_mode == LayoutMode::SoAPC) {}
    else
        for (auto const& p : arr)
            out << p.copy();
    return out;
}




template<auto opts, typename internals>
void empty(ParticleArray<opts, internals>& array)
{
    array.clear();
}


template<auto opts, typename internals>
void swap(ParticleArray<opts, internals>& array1, ParticleArray<opts, internals>& array2)
{
    array1.swap(array2);
}

template<typename P0, typename P1>
EqualityReport particle_compare(P0 const& p0, P1 const& p1, std::size_t const i = 0,
                                double const atol = 1e-15)
{
    std::string idx = std::to_string(i);
    if (p0.iCell() != p1.iCell())
        return EqualityReport{false, "icell mismatch at index: " + idx, i};

    if (!float_equals(p0.v(), p1.v(), atol))
        return EqualityReport{false, "v mismatch at index: " + idx, i};
    if (!float_equals(p0.delta(), p1.delta(), atol))
        return EqualityReport{false, "delta mismatch at index: " + idx, i};


    return EqualityReport{true};
}


template<typename P0, typename P1>
EqualityReport particles_equals(P0 const& ref, P1 const& cmp, double const atol = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ref.size()) + " vs "
                                         + std::to_string(cmp.size())};

    auto rit      = ref.begin();
    auto cit      = cmp.begin();
    std::size_t i = 0;

    for (; rit != ref.end(); ++rit, ++cit, ++i)
        if (auto const eq = particle_compare(*rit, *cit, i, atol); !eq)
            return eq;

    return EqualityReport{true};
}

template<auto O, typename I0, typename I1>
EqualityReport operator==(ParticleArray<O, I0> const& p0, ParticleArray<O, I1> const& p1)
{
    auto report = particles_equals(p0, p1);
    if (!report)
    {
        PHARE_LOG_LINE_STR(p0[report.idx].copy());
        PHARE_LOG_LINE_STR(p1[report.idx].copy());
    }
    return report;
}

template<auto O, typename I0>
EqualityReport operator==(ParticleArray<O, I0> const& p0, std::vector<Particle<O>> const& p1)
{
    return particles_equals(p0, p1);
}



template<typename ParticleArray_t>
auto constexpr base_layout_type()
{
    using enum LayoutMode;
    if constexpr (any_in(ParticleArray_t::layout_mode, AoS, AoSMapped, AoSPC, AoSTS))
        return AoS;
    return SoA;
}


template<typename ParticleArray_t>
auto count_particles_per_cell(ParticleArray_t const& ps)
{
    std::unordered_map<std::string, std::size_t> count_per_cell;

    auto inc_cell = [&](auto const& cell) {
        auto const p = Point{cell}.str();
        if (count_per_cell.count(p) == 0)
            count_per_cell.emplace(p, 0);
        ++count_per_cell[p];
    };

    if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS))
        for (auto const& tile : ps())
            for (auto const& p : tile())
                inc_cell(p.iCell());
    else
        for (auto const& p : ps)
            inc_cell(p.iCell());

    return count_per_cell;
}

template<typename ParticleArray_t>
auto print_particles_per_cell(ParticleArray_t const& ps)
{
    for (auto [k, v] : count_particles_per_cell(ps))
    {
        PHARE_LOG_LINE_SS(k << " " << v);
    }
}


template<auto O, typename I0>
void per_particle(ParticleArray<O, I0> const& particles, auto const fn)
{
    using enum LayoutMode;
    if constexpr (any_in(O.layout_mode, AoSTS))
        for (auto& tile : particles())
            for (auto& p : tile())
                fn(p);
    else
        for (auto& p : particles)
            fn(p);
}

template<auto O, typename I0>
void check_particles(ParticleArray<O, I0> const& particles, bool const print = false)
{
    using enum LayoutMode;

    PHARE_DEBUG_DO({
        if constexpr (any_in(O.layout_mode, AoSTS))
            particles.check();
    })
}


template<auto O, typename I0>
void check_particles_views(ParticleArray<O, I0>& particles, bool const print = false)
{
    using enum LayoutMode;

    PHARE_DEBUG_DO({ check_particles(*particles); })
}

template<auto O, typename I0>
void particle_array_domain_is_valid(ParticleArray<O, I0> const& particles, auto const& domain_box)
{
    std::size_t in_domain_box = 0, not_in_domain_box = 0;

    if constexpr (any_in(O.layout_mode, AoSTS))
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto const& tile = particles()[tidx];
            for (std::size_t i = 0; i < tile().size(); ++i)
                if (not isIn(tile()[i], tile))
                    throw std::runtime_error("Invalid location");
        }
    }

    per_particle(particles, [&](auto const& p) {
        if (isIn(p, domain_box))
            ++in_domain_box;
        else
            ++not_in_domain_box;
    });

    if (not(not_in_domain_box == 0 and in_domain_box == particles.size()))
        throw std::runtime_error("Invalid particles");
}

template<auto O, typename I0>
void particle_array_ghost_is_valid(ParticleArray<O, I0> const& particles, auto const& domain_box,
                                   auto const& ghost_box)
{
    std::size_t in_ghost_layer = 0, not_in_ghost_layer = 0, outside_gb = 0;

    per_particle(particles, [&](auto const& p) {
        auto const ingb = isIn(p, ghost_box);
        auto const indb = isIn(p, domain_box);

        if (!ingb)
            ++outside_gb;
        else if (ingb and not indb)
            ++in_ghost_layer;
        else
            ++not_in_ghost_layer;
    });

    if constexpr (any_in(O.layout_mode, AoSTS))
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto const& tile = particles()[tidx];
            for (std::size_t i = 0; i < tile().size(); ++i)
                if (not isIn(tile()[i], tile))
                    throw std::runtime_error("Invalid location");
        }
    }

    if (not(outside_gb == 0 and not_in_ghost_layer == 0 and in_ghost_layer == particles.size()))
        throw std::runtime_error("Invalid particles");
}


} // namespace PHARE::core


#endif
