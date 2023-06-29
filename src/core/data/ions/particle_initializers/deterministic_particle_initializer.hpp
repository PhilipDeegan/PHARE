#ifndef PHARE_CORE_DETERMINISTIC_PARTICLE_INITIALIZER_HPP
#define PHARE_CORE_DETERMINISTIC_PARTICLE_INITIALIZER_HPP

#include <memory>
#include <cassert>
#include <functional>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE::core
{
template<typename ParticleArray, typename GridLayout>
class DeterministicParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
{
    using Particle = typename ParticleArray::value_type;

public:
    static constexpr auto dimension = GridLayout::dimension;
    using InputFunction             = initializer::InitFunction<dimension>;

    DeterministicParticleInitializer(InputFunction const& density,
                                     std::array<InputFunction, 3> const& bulkVelocity,
                                     std::array<InputFunction, 3> const& thermalVelocity,
                                     double const& particleCharge,
                                     std::uint32_t const& nbrParticlesPerCell)
        : density_{density}
        , bulkVelocity_{bulkVelocity}
        , thermalVelocity_{thermalVelocity}
        , particleCharge_{particleCharge}
        , nbrParticlePerCell_{nbrParticlesPerCell}
    {
    }
    virtual ~DeterministicParticleInitializer() = default;
    void loadParticles(ParticleArray& particles, GridLayout const& layout) const override;

private:
    InputFunction const density_;
    std::array<InputFunction, 3> const bulkVelocity_;
    std::array<InputFunction, 3> const thermalVelocity_;

    double const particleCharge_;
    std::uint32_t const nbrParticlePerCell_;
};


class DeterministicInitFunctions
{
public:
    template<typename Function, typename FunctionArray, typename... Coords>
    DeterministicInitFunctions(Function& density, FunctionArray& bulkVelocity,
                               FunctionArray& thermalVelocity, Coords const&... coords)
        : _n{density(coords...)}
    {
        static_assert(sizeof...(coords) <= 3, "can only provide up to 3 coordinates");
        for (std::uint32_t i = 0; i < 3; i++)
        {
            _V[i]   = bulkVelocity[i](coords...);
            _Vth[i] = thermalVelocity[i](coords...);
        }
    }

    std::array<double const*, 3> B() const { return ptrs(_B); }

    auto operator()() const { return std::make_tuple(_n->data(), ptrs(_V), ptrs(_Vth)); }

private:
    std::array<double const*, 3>
    ptrs(std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> const& v) const
    {
        return {v[0]->data(), v[1]->data(), v[2]->data()};
    }

    std::shared_ptr<PHARE::core::Span<double>> const _n;
    std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> _B, _V, _Vth;
};

template<std::size_t ndim>
std::vector<std::array<double, ndim>> velocity_directions()
{
    if constexpr (ndim == 1) // left / right
        return {{-1}, {1}};
    if constexpr (ndim == 2) // eight directions including diagonals
        return {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
    if constexpr (ndim == 3) // 24 directions?
        return {{-1, -1, -1}, {-1, 0, -1}, {-1, 1, -1}, {0, -1, -1}, {0, 1, -1},  {1, -1, -1},
                {1, 0, -1},   {1, 1, -1},  {-1, -1, 0}, {-1, 0, 0},  {-1, 1, 0},  {0, -1, 0},
                {0, 1, 0},    {1, -1, 0},  {1, 0, 0},   {1, 1, 0},   {-1, -1, 1}, {-1, 0, 1},
                {-1, 1, 1},   {0, -1, 1},  {0, 1, 1},   {1, -1, 1},  {1, 0, 1},   {1, 1, 1}};
    throw std::runtime_error("FAIL");
}
template<std::size_t ndim>
void deterministicVelocity(std::array<double, 3> const& V, std::array<double, 3> const& /*Vth*/,
                           std::array<double, ndim> const& direction,
                           std::array<double, 3>& partVelocity)
{
    partVelocity[0] = V[0] * direction[0];
    if constexpr (ndim > 1)
        partVelocity[1] = V[0] * direction[1];
    if constexpr (ndim > 2)
        partVelocity[2] = V[0] * direction[2];
}

template<typename ParticleArray, typename GridLayout>
void DeterministicParticleInitializer<ParticleArray, GridLayout>::loadParticles(
    ParticleArray& particles, GridLayout const& layout) const
{
    auto point = [](std::size_t i, auto const& indices) -> core::Point<std::uint32_t, dimension> {
        if constexpr (dimension == 1)
            return {std::get<0>(indices[i])};
        if constexpr (dimension == 2)
            return {std::get<0>(indices[i]), std::get<1>(indices[i])};
        if constexpr (dimension == 3)
            return {std::get<0>(indices[i]), std::get<1>(indices[i]), std::get<2>(indices[i])};
    };
    auto deltas        = ConstArray<double, dimension>(.5); // all start in the middle
    auto ndCellIndices = layout.physicalStartToEndIndices(QtyCentering::primal);
    auto cellCoords    = layout.indexesToCoordVectors(
           ndCellIndices, QtyCentering::primal, [](auto const& gridLayout, auto const&... indexes) {
            return gridLayout.cellCenteredCoordinates(indexes...);
        });
    auto const fns         = std::make_from_tuple<DeterministicInitFunctions>(std::tuple_cat(
                std::forward_as_tuple(density_, bulkVelocity_, thermalVelocity_), cellCoords));
    auto const [n, V, Vth] = fns();
    auto directions        = velocity_directions<dimension>();
    for (std::size_t flatCellIdx = 0; flatCellIdx < ndCellIndices.size(); flatCellIdx++)
    {
        auto const cellWeight   = n[flatCellIdx] / nbrParticlePerCell_;
        auto const AMRCellIndex = layout.localToAMR(point(flatCellIdx, ndCellIndices));

        for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
        {
            std::array<double, 3> particleVelocity{0, 0, 0};
            deterministicVelocity(
                {V[0][flatCellIdx], V[1][flatCellIdx], V[2][flatCellIdx]},
                {Vth[0][flatCellIdx], Vth[1][flatCellIdx], Vth[2][flatCellIdx]}, //
                directions[ipart % directions.size()], particleVelocity);

            particles.emplace_back(Particle{cellWeight, particleCharge_,
                                            AMRCellIndex.template toArray<int>(), deltas,
                                            particleVelocity});
        }
    }
}

} // namespace PHARE::core


#endif /* PHARE_CORE_DETERMINISTIC_PARTICLE_INITIALIZER_HPP  */
