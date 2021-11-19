
#include "bench/real/bench_harris.h"

namespace PHARE::amr::bench
{
template<std::size_t dim, std::size_t interp, std::size_t nbRefinePart>
struct HarrisBench : public benchmark::Fixture
{
    using Simulator_t = PHARE::Simulator<dim, interp, nbRefinePart>;

public:
    void SetUp(::benchmark::State const& state) override
    {
        if (!hierarchy)
        {
            auto& dict = PHARE::real::bench::harris::createDict();
            hierarchy  = PHARE::amr::Hierarchy::make();
            simulator  = std::make_unique<Simulator_t>(dict, hierarchy);
            simulator->initialize();
            simulator->dump(simulator->currentTime(), simulator->timeStep());
            PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        }
    }

    void run(::benchmark::State&);
    void TearDown(::benchmark::State const& /*state*/) override
    {
        simulator.reset();
        hierarchy.reset();
        PHARE::SamraiLifeCycle::reset();
    }

    std::shared_ptr<Hierarchy> hierarchy;
    std::unique_ptr<Simulator_t> simulator;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void HarrisBench<dim, interp, op>::run(::benchmark::State& state)
{
    for (auto _ : state)
        while (simulator->currentTime() < simulator->endTime())
        {
            simulator->advance(simulator->timeStep());
            simulator->dump(simulator->currentTime(), simulator->timeStep());
        }
}

BENCHMARK_TEMPLATE_DEFINE_F(HarrisBench, _2_1_advance, 2, 1, 4)(benchmark::State& state)
{
    run(state);
}
BENCHMARK_REGISTER_F(HarrisBench, _2_1_advance)->Unit(benchmark::kNanosecond);
// ->ArgsProduct({THREADS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
