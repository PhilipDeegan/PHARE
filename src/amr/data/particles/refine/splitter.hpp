

#ifndef PHARE_SPLITTER_HPP
#define PHARE_SPLITTER_HPP

#include <array>
#include <cmath>
#include <tuple>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "amr/amr_constants.hpp"

namespace PHARE::amr
{
// see meta_utilities.hpp for list of declared permutations.

template<typename _dimension, typename _nbRefinedPart>
struct SplitPattern
{
    constexpr SplitPattern(float weight)
        : weight_{weight}
    {
    }

    float weight_;
    std::array<core::Point<float, _dimension{}()>, _nbRefinedPart{}()> deltas_{};
};


template<typename... Patterns>
class PatternDispatcher
{
public:
    constexpr PatternDispatcher(Patterns&&... _patterns)
        : patterns{_patterns...}
        , nbRefinedParts{(_patterns.deltas_.size() + ...)}
    {
    }

    template<typename FineiCellDelta, typename ParticleIt, typename Particles>
    inline void operator()(FineiCellDelta const& fine, ParticleIt const& coarseParticle,
                           Particles& refinedParticles, size_t idx = 0) const
    {
        dispatch(fine, coarseParticle, refinedParticles, idx);
    }

    std::tuple<Patterns...> patterns{};
    size_t nbRefinedParts{0};

private:
    template<typename FineiCellDelta, typename ParticleIt, typename Particles>
    void dispatch(FineiCellDelta const& fine, ParticleIt const& particle, //
                  Particles& particles, size_t idx) const
    {
        using Weight_t = std::decay_t<decltype(particle.weight())>;
        using Delta_t  = std::decay_t<decltype(particle.delta()[0])>;

        constexpr auto dimension = ParticleIt::dimension;
        constexpr auto refRatio  = PHARE::amr::refinementRatio;
        constexpr std::array power{refRatio, refRatio * refRatio, refRatio * refRatio * refRatio};

        assert(particles.size() >= idx + nbRefinedParts);

        auto split = [&](auto& fineParticle, auto& pattern, auto const& rpIndex) {
            fineParticle.weight()
                = particle.weight() * static_cast<Weight_t>(pattern.weight_) * power[dimension - 1];
            fineParticle.charge() = particle.charge();
            fineParticle.iCell()  = fine.iCell();
            fineParticle.delta()  = fine.delta();
            fineParticle.v()      = particle.v();

            for (size_t iDim = 0; iDim < dimension; iDim++)
            {
                fineParticle.delta()[iDim] += static_cast<Delta_t>(pattern.deltas_[rpIndex][iDim]);
                Delta_t integra = std::floor(fineParticle.delta()[iDim]);
                fineParticle.delta()[iDim] -= integra;
                fineParticle.iCell()[iDim] += static_cast<int32_t>(integra);
            }
        };

        core::apply(patterns, [&](auto const& pattern) {
            for (size_t rpIndex = 0; rpIndex < pattern.deltas_.size(); rpIndex++)
            {
                auto fineParticle = particles.begin() + idx++;

                if constexpr (Particles::layout_mode == core::LayoutMode::SoA)
                    split(fineParticle, pattern, rpIndex);
                else
                    split(*fineParticle, pattern, rpIndex);
            }
        });
    }
};


template<typename _dimension, typename _interp_order, typename _nbRefinedPart>
struct ASplitter
{
    static constexpr auto dimension     = _dimension{}();
    static constexpr auto interp_order  = _interp_order{}();
    static constexpr auto nbRefinedPart = _nbRefinedPart{}();

    static constexpr int maxCellDistanceFromSplit()
    {
        constexpr auto particleSize = interp_order + 1;
        return std::ceil(particleSize * 0.5);
    }

    constexpr ASplitter() {}
};

template<typename dimension, typename interp_order, typename nbRefinedPart>
class Splitter : public ASplitter<dimension, interp_order, nbRefinedPart>
{
    Splitter() = delete; // Unspecialized template class, never to be instantiated
};

template<typename dim>
struct BlackDispatcher : SplitPattern<dim, core::RefinedParticlesConst<1>>
{
    using Super = SplitPattern<dim, core::RefinedParticlesConst<1>>;
    constexpr BlackDispatcher(float const weight)
        : Super{weight}
    {
    }
};

template<typename dim>
struct PurpleDispatcher
{
};
template<typename dim>
struct BrownDispatcher
{
};
template<typename dim>
struct PinkDispatcher
{
};

} // namespace PHARE::amr

#endif // endif PHARE_SPLITTER_HPP
