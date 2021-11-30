
#include "mkn/kul/dbg.hpp"
#include "mkn/kul/log.hpp"

#include "bench/amr/solver/solver_ppc_bench.h"

#include <atomic>
#include <thread>

#include "core/numerics/ohm/ohm.h"
#include "core/numerics/ampere/ampere.h"
#include "core/numerics/faraday/faraday.h"
#include "core/numerics/ion_updater/ion_range_updater.h"


namespace PHARE::amr::bench
{
void abort_if(bool b)
{
    if (b)
        std::abort();
}

std::size_t static constexpr thread_start = _PHARE_BENCH_THREAD_START_;
std::size_t static constexpr thread_end   = _PHARE_BENCH_THREADS_END_;

static auto threads()
{
    std::vector<std::int64_t> vec;
    for (std::size_t i = thread_start; i <= thread_end; ++i)
        vec.emplace_back(i);
    return vec;
}
static std::vector<std::int64_t> THREADS{/*threads()*/ 10};

template<std::size_t dim, std::size_t interp, std::size_t op = 1>
struct SolveFixture : public benchmark::Fixture
{
    using PHARE_Types = PHARE::core::PHARE_Types<dim, interp>;

    using ParticleArray_t = typename PHARE_Types::ParticleArray_t;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using Ions_t          = typename PHARE_Types::Ions_t;
    using PatchState      = typename PHARE::core::bench::HybridPatch<GridLayout_t>::State;

    using Electromag_t = PHARE::core::bench::Electromag<GridLayout_t>;
    using IonUpdater   = PHARE::core::IonRangeUpdater<Ions_t, Electromag_t, GridLayout_t>;

public:
    void SetUp(::benchmark::State const& state) override
    {
        constexpr std::uint32_t cells = 1e3;
        constexpr std::uint32_t parts = 1e7;
        assert(state.range(0) > 0);
        std::uint32_t n_threads = static_cast<std::uint32_t>(state.range(0));
        assert(n_threads);
        for (std::size_t i = 0; i < n_threads; ++i)
            patches.emplace_back(PatchState::make_unique(
                PHARE::core::ConstArray<std::uint32_t, dim>(0),
                PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), parts));
        assert(patches.size());
    }

    void TearDown(::benchmark::State const& /*state*/) override { patches.clear(); }

    void solve(::benchmark::State&);

    std::vector<std::unique_ptr<PatchState>> patches;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void SolveFixture<dim, interp, op>::solve(::benchmark::State& state)
{
    assert(patches.size());
    auto n_threads_ = state.range(0);
    abort_if(n_threads_ < 1);
    std::uint16_t n_threads = n_threads_;

    auto units = PHARE::amr::bench::solver_update_dao_per_thread<PatchState>(patches, n_threads);

    PHARE::core::RangeSynchrotron<ParticleArray_t> synchrotron{n_threads}; // syncs on destruct

    auto per_thread = [&](std::uint16_t idx) {
        IonUpdater{synchrotron, idx, "modified_boris"}.updatePopulations(
            units[idx], 1e-5, PHARE::core::UpdaterMode::domain_only);
    };

    for (auto _ : state)
    {
        auto threads = PHARE::core::generate(
            [&](auto i) { return std::thread{[&, i]() { per_thread(i); }}; }, n_threads - 1);
        per_thread(n_threads - 1);

        for (auto& thread : threads)
            if (thread.joinable())
                thread.join();
    }
}


BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_push, 2, 1, 2)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_push)->Unit(benchmark::kNanosecond)->ArgsProduct({THREADS});
} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
