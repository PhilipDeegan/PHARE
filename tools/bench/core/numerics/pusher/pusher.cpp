
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tools/bench/core/bench.hpp"

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    using namespace PHARE::core;
    constexpr double mass           = 1;
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;
    auto constexpr no_op            = [](auto& particleRange) { return particleRange; };

    using PHARE_Types     = PHARE_Types<dim, interp>;
    using GridLayout_t    = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator    = Interpolator<dim, interp>;
    using ParticleArray_t = typename PHARE_Types::ParticleArray_t;
    using Electromag_t    = UsableElectromag<GridLayout_t, ParticleArray_t::alloc_mode>;
    using BorisPusher_t = BorisPusher<dim, IndexRange<ParticleArray_t>, Electromag_t, Interpolator,
                                      BoundaryCondition<dim, interp>, GridLayout_t>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    ParticleArray_t domainParticles{layout.AMRBox(), n_parts, particle<dim>()};
    disperse(domainParticles, layout.AMRBox(), 13337);
    ParticleArray_t tmpDomain = domainParticles;
    auto rangeIn              = makeIndexRange(domainParticles);
    auto rangeOut             = makeIndexRange(tmpDomain);
    BorisPusher_t pusher{layout.meshSize(), .001};
    Interpolator interpolator;

    while (state.KeepRunningBatch(1))
        pusher.move(rangeIn, rangeOut, em, mass, interpolator, layout, no_op, no_op);
}

BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
