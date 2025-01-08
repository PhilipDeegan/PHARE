#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/particles/particle.hpp"

#include <array>
#include <cstdint>

namespace PHARE::core
{

enum class LayoutMode : std::uint16_t {
    AoS = 0,
    AoSMapped, // 1
    AoSPC,     // 2
    AoSTS,     // 3
    SoA,       // 4
    SoAVX,     // 5
    SoATS,     // 6
    SoAVXTS,   // 7
    SoAPC      // 8
};

enum class StorageMode : std::uint16_t { ARRAY = 0, VECTOR, SPAN };

enum class ParticleType : std::uint16_t { Domain = 0, Ghost, PatchGhost, LevelGhost, All };


template<std::size_t dim>
struct ParticleDefaults
{
    using Particle_t = Particle<dim>;
};


template<typename R = std::uint32_t, std::size_t dim>
auto as_local_cell(std::array<int, dim> const& lower,
                   std::array<int, dim> const& icell) _PHARE_ALL_FN_
{
    return array_minus<R>(icell, lower);
}
template<std::size_t dim>
auto as_local_cell(Box<int, dim> const& box, std::array<int, dim> const& icell) _PHARE_ALL_FN_
{
    return as_local_cell(box.lower.toArray(), icell);
}
template<std::size_t dim>
auto as_local_cell(Box<int, dim> const& box, Point<int, dim> const& icell) _PHARE_ALL_FN_
{
    return as_local_cell(box.lower.toArray(), icell.toArray());
}

template<typename I0, typename I1> // support const vs non-const iterators
auto it_dist(I0&& i0, I1&& i1)
{
    auto const& begin = i0;
    auto const& pos   = i1;
    auto d            = std::distance(begin, pos);
    PHARE_ASSERT(d < 1e18); // should never happen // eg 2635249153387078728
    return d;
}



template<typename Box_t, typename RValue = std::uint32_t>
class LocalisedCellFlattener
{
public:
    static constexpr std::size_t dim = Box_t::dimension;

    LocalisedCellFlattener(Box_t const& b)
        : box{b}
        , shape{box.shape()}
    {
    }

    template<typename int_t>
    RValue operator()(std::array<int_t, dim> const icell) const
    {
        for (std::size_t i = 0; i < dim; ++i)
            icell[i] -= box.lower[i];
        if constexpr (dim == 2)
            return icell[1] + icell[0] * shape[1];
        if constexpr (dim == 3)
            return icell[2] + icell[1] * shape[2] + icell[0] * shape[1] * shape[2];
        return icell[0];
    }
    template<typename Particle>
    RValue operator()(Particle const& particle) const
    {
        return (*this)(particle.iCell);
    }

    Box_t const& box;

private:
    Point<int, dim> const& shape;
};



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP*/
