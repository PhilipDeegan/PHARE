#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTING_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTING_HPP

#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/sorting.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "particle.hpp"
#include "particle_array.hpp"

#include <iostream>
#include <cstddef>
#include <vector>

namespace PHARE::core
{

template<typename box_, typename T>
struct BoxSpan : Span<T>
{
    using box_t = box_;

    box_t box;
};

template<typename ParticleArray_t>
struct ParticleCountSorting
{
    auto constexpr static dim = ParticleArray_t::dimension;
    using box_t               = typename ParticleArray_t::box_t;

    void operator()()
    {
        sorter.setup(particles.size(), particles.box());
        sorter.sort(particles, flattener);
        cell_converter(particles.box(), flattener);
    }

    ParticleArray_t& particles;
    CountingSort<ParticleArray_t, box_t::dimension>& sorter;
    CountingSortCellConverter<ParticleArray_t> cell_converter{sorter, particles};
    LocalisedCellFlattener<box_t> flattener{particles.box()};
};

template<typename ParticleArray_t>
struct ParticleCountRangeFinder
{
    auto constexpr static dim = ParticleArray_t::dimension;
    using box_t               = typename ParticleArray_t::box_t;
    using Particle_t          = typename ParticleArray_t::value_type;
    using Span_t              = BoxSpan<box_t, Particle_t const>;
    using Ranges_t            = std::vector<Span_t>;

    auto ranges(box_t const& box) const;
    auto copy(box_t const& box) const;
    auto copy_into(box_t const& box, ParticleArray_t& into) const;

    ParticleArray_t const& particles;
};



template<typename ParticleArray_t>
struct ParticleRangeEraser
{
    auto constexpr static dim = ParticleArray_t::dimension;
    using box_t               = typename ParticleArray_t::box_t;
    using Particle_t          = typename ParticleArray_t::value_type;

    using Span_t   = BoxSpan<box_t, Particle_t const>;
    using Ranges_t = std::vector<Span_t>;

    void erase(box_t const& box);
    void erase_outside(box_t const& box);

    ParticleCountSorting<ParticleArray_t>& sorter;
    std::vector<decltype(Particle_t{}.iCell)> iCells{};
    ParticleCountRangeFinder<ParticleArray_t> finder{sorter.particles};
};



template<typename ParticleArray_t>
auto ParticleCountRangeFinder<ParticleArray_t>::ranges(box_t const& box) const
{
    // PHARE_LOG_LINE_STR(box.lower << " " << box.upper);

    auto begin = particles.vector().data();

    auto lo = box.lower.toArray();
    auto up = box.upper.toArray();
    Ranges_t ranges;

    auto add = [&](auto const&... args) {
        std::array<std::int32_t, dim> l{args..., lo[dim - 1]};
        box_t box{l, {args..., up[dim - 1]}};
        ranges.emplace_back(Span_t{begin + particles.ppc_offsets(l),
                                   sum_from(
                                       [&](/*auto=FAIL?*/ Point<std::int32_t, dim> const& cell) {
                                           return particles.ppc(cell);
                                       },
                                       box),
                                   box});
    };

    if constexpr (dim == 1)
        add();

    else if constexpr (dim == 2)
        for (auto ix = lo[0]; ix <= up[0]; ++ix)
            add(ix);

    else if constexpr (dim == 3)
        for (auto ix = lo[0]; ix <= up[0]; ++ix)
            for (auto iy = lo[1]; iy <= up[1]; ++iy)
                add(ix, iy);

    return ranges;
}

template<typename ParticleArray_t>
auto ParticleCountRangeFinder<ParticleArray_t>::copy_into(box_t const& box,
                                                          ParticleArray_t& copy) const
{
    auto ranges = this->ranges(box);
    copy.reserve(copy.size() + sum_from([](auto const& range) { return range.size(); }, ranges));
    for (auto const& range : ranges)
        copy.insert(copy.end(), range.begin(), range.end());
    return copy;
}

template<typename ParticleArray_t>
auto ParticleCountRangeFinder<ParticleArray_t>::copy(box_t const& box) const
{
    ParticleArray_t copy{box};
    copy_into(box, copy);
    return copy;
}


template<typename ParticleArray_t>
void ParticleRangeEraser<ParticleArray_t>::erase(box_t const& erase_box)

{
    using PartIterator = typename ParticleArray_t::const_iterator;

    for (auto& range : reverse(finder.ranges(erase_box)))
        sorter.particles.erase(PartIterator{range.begin()}, PartIterator{range.end()});
}

template<typename ParticleArray_t>
void ParticleRangeEraser<ParticleArray_t>::erase_outside(box_t const& keep_box)

{
    iCells.clear(); // zero size
    for (auto outside_keep_box : sorter.particles.box().remove(keep_box))
        for (auto const& range : finder.ranges(outside_keep_box))
            for (auto const& cell : range.box)
                iCells.emplace_back(cell.toArray());

    std::sort(iCells.begin(), iCells.end(), [&](auto const& a, auto const& b) {
        return sorter.flattener(a) > sorter.flattener(b);
    }); // last to first

    for (auto& cell : iCells)
    {
        auto begin = sorter.particles.begin() + sorter.particles.ppc_offsets(cell);
        auto end   = begin + sorter.particles.ppc(cell);
        sorter.particles.erase(begin, end);
    }
}



} // namespace PHARE::core
#endif
