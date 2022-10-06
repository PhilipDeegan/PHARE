

#include "bench/core/bench.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/faraday/faraday_avx.hpp"

constexpr std::size_t dim     = 3;
constexpr std::size_t interp  = 1;
constexpr std::uint32_t cells = 500;
constexpr double dt           = .001;

using PHARE_Types  = PHARE::core::PHARE_Types<dim, interp>;
using GridLayout_t = typename PHARE_Types::GridLayout_t;
using Faradat_y    = PHARE::core::Faraday<GridLayout_t>;
using PHARE::core::Component;
using PHARE::core::Direction;

GridLayout_t getLayout()
{
    return {PHARE::core::ConstArray<double, dim>(1.0 / cells),
            PHARE::core::ConstArray<std::uint32_t, dim>(cells),
            PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}};
}

void deriv(benchmark::State& state)
{
    auto layout = getLayout();
    PHARE::core::bench::Electromag EM{layout}, EMNew{layout};
    PHARE::core::Faraday<GridLayout_t> faraday;
    auto __ = PHARE::core::SetLayout(&layout, faraday);

    while (state.KeepRunning())
        faraday(EM.B, EM.E, EMNew.B, dt);
}


void derivX(benchmark::State& state)
{
    auto layout = getLayout();
    PHARE::core::FaradayHelper::I().set(layout); // alloc before bench
    PHARE::core::bench::Electromag EM{layout}, EMNew{layout};
    PHARE::core::FaradayVX<GridLayout_t> faraday;
    auto __ = PHARE::core::SetLayout(&layout, faraday);

    while (state.KeepRunning())
        faraday(EM.B, EM.E, EMNew.B, dt);
}



BENCHMARK(derivX)->Unit(benchmark::kMicrosecond);
BENCHMARK(deriv)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
