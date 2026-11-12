// MultiBoris push / deposit benchmark over a level of patches: only MultiBoris::move is timed
//  (COPY: push copies + deposit, REF: in place push incl. tiled layouts' cell tracking), no
//  sorting, partitioning or erasing of particles leaving patches: REF iterations run the updater's
//  post_move_sync before each move, untimed, as moves are not valid back to back without it
//
// the level: a scenario file (PHARE_BENCH_SCENARIO), one patch per line
//    level ppc lower_0 .. lower_dim-1 upper_0 .. upper_dim-1
//  "#" starts a comment, "# dl x [y [z]]" sets the mesh size (default 0.1)
//  or without a scenario: PHARE_PATCHES cubes of PHARE_CELLS cells, PHARE_PPC particles per cell
//
// env:
//  PHARE_BENCH_SCENARIO  scenario file, see above
//  PHARE_BENCH_LEVEL     the level of the scenario benchmarked (default 0)
//  PHARE_BENCH_SCALE     boxes shrunk by this factor per direction, particles per cell kept (1)
//  PHARE_BENCH_MIN_CELLS smallest extent per direction of a shrunk box (default 1)
//  PHARE_BENCH_V         velocities in [-v, v], default 0: particles stay put so REF repeats as is
//  PHARE_SEED            particle generation, per patch seed + patch index (default 1337)
//  PHARE_PATCHES, PHARE_CELLS, PHARE_PPC: synthetic level, see above (10 patches of 20^dim, 100)
//  PHARE_BENCH_POOLS, PHARE_BENCH_THREADS: CPU thread pools x threads per pool (default 1 x 1)
//  GPU: one mkn.gpu ThreadedStreamLauncher host thread, a stream per patch, its per patch
//   function times dumped to $PHARE_ASYNC_TIMES (see launcher_times), move_ns / sync_ns per patch
//
// phase split: phlop scopes in MultiBoris (-DPHARE_LOG_LEVEL=3, env PHARE_SCOPE_TIMING=1), and
//  NVTX ranges with -DPHLOP_NVTX for nsys

#include "core/logger.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/data/particles/particle_array_appender.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"
#include "tests/amr/resources_manager/test_resources_manager_fixtures.hpp"

#include "benchmark/benchmark.h"

#include <array>
#include <string>
#include <vector>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <stdexcept>

namespace PHARE::core::bench
{
auto static const scenario  = get_env_as("PHARE_BENCH_SCENARIO", std::string{});
auto static const level     = get_env_as("PHARE_BENCH_LEVEL", int{0});
auto static const scale     = get_env_as("PHARE_BENCH_SCALE", int{1});
auto static const min_cells = get_env_as("PHARE_BENCH_MIN_CELLS", int{1});
auto static const v_max     = get_env_as("PHARE_BENCH_V", double{0});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1337});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{10});
auto static const cells     = get_env_as("PHARE_CELLS", int{20});
auto static const synth_ppc = get_env_as("PHARE_PPC", std::size_t{100});
double static constexpr dt  = .001;

// CPU MultiBoris dispatch: pools x threads per pool, set before the pool's first use
bool static const thread_pool_setup = []() {
    ThreadPool::n_pools          = get_env_as("PHARE_BENCH_POOLS", std::size_t{1});
    ThreadPool::threads_per_pool = get_env_as("PHARE_BENCH_THREADS", std::size_t{1});
    return true;
}();


template<std::size_t dim>
struct ScenarioPatch
{
    Box<int, dim> box;
    std::size_t ppc;
};

template<std::size_t dim>
struct Scenario
{
    std::array<double, dim> dl = ConstArray<double, dim>(.1);
    std::vector<ScenarioPatch<dim>> patches;
};


// box shrunk by scale per direction around its lower corner, at least min_cells
template<std::size_t dim>
Box<int, dim> scaled(Box<int, dim> const& box)
{
    auto out = box;
    for (std::size_t d = 0; d < dim; ++d)
    {
        int const extent = box.upper[d] - box.lower[d] + 1;
        int const shrunk = std::max(min_cells, (extent + scale - 1) / scale);
        out.upper[d]     = box.lower[d] + shrunk - 1;
    }
    return out;
}


template<std::size_t dim>
Scenario<dim> load_scenario()
{
    Scenario<dim> scen;

    if (scenario.empty()) // synthetic: cubes next to each other in x
    {
        for (std::size_t i = 0; i < n_patches; ++i)
        {
            Box<int, dim> box;
            for (std::size_t d = 0; d < dim; ++d)
            {
                box.lower[d] = d == 0 ? static_cast<int>(i) * cells : 0;
                box.upper[d] = box.lower[d] + cells - 1;
            }
            scen.patches.push_back({scaled(box), synth_ppc});
        }
        return scen;
    }

    std::ifstream file{scenario};
    if (!file)
        throw std::runtime_error("PHARE_BENCH_SCENARIO: cannot open " + scenario);

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream in{line};
        std::string first;
        if (!(in >> first))
            continue;
        if (first == "#")
        {
            std::string key;
            if (in >> key and key == "dl")
                for (std::size_t d = 0; d < dim; ++d)
                    in >> scen.dl[d];
            continue;
        }
        if (first[0] == '#')
            continue;

        int const lvl = std::stoi(first);
        std::size_t ppc;
        Box<int, dim> box;
        in >> ppc;
        for (std::size_t d = 0; d < dim; ++d)
            in >> box.lower[d];
        for (std::size_t d = 0; d < dim; ++d)
            in >> box.upper[d];
        if (!in)
            throw std::runtime_error("PHARE_BENCH_SCENARIO: bad line for dim " + std::to_string(dim)
                                     + ": " + line);
        if (lvl == level)
            scen.patches.push_back({scaled(box), ppc});
    }

    if (scen.patches.empty())
        throw std::runtime_error("PHARE_BENCH_SCENARIO: no patch on level "
                                 + std::to_string(level));
    return scen;
}


template<auto opts>
auto make_layout(auto const& dl, auto const& box)
{
    using GridLayout_t = PHARE_Types<opts>::Hybrid::GridLayout_t;
    auto constexpr dim = opts.dimension;
    Point<double, dim> origin;
    for (std::size_t d = 0; d < dim; ++d)
        origin[d] = box.lower[d] * dl[d];
    return GridLayout_t{dl, box.shape().template toArray<std::uint32_t>(), origin, box};
}


// one patch's particles: generated as AoSMapped (the test fixtures), appended to the layout
template<auto opts>
auto make_ions(auto const& layout, std::size_t const ppc, std::size_t const patch_idx)
{
    auto constexpr static ref_opts = SimOpts{.dimension = opts.dimension};
    auto const patch_seed          = static_cast<int>(seed + patch_idx);

    UsableIons_t<ref_opts> ref{make_layout<ref_opts>(layout.meshSize(), layout.AMRBox()),
                               "protons"};
    {
        auto& particles = ref.populations[0].particles;
        add_particles(particles.domain_particles, layout.AMRBox(), ppc);
        add_ghost_particles(particles.level_ghost_particles, layout.AMRBox(), ppc, 1);
        for (auto* ps : {&particles.domain_particles, &particles.level_ghost_particles})
        {
            delta_disperse(*ps, patch_seed);
            vary_velocity(*ps, -v_max, v_max, patch_seed);
        }
    }

    if constexpr (opts.layout_mode == LayoutMode::AoSMapped) // generated as such already
        return ref;
    else
    {
        UsableIons_t<opts> ions{layout, "protons"};
        auto& src = ref.populations[0].particles;
        auto& dst = ions.populations[0].particles;
        reserve(dst.domain_particles, layout.AMRBox(), ppc);
        append_particles<ParticleType::Domain>(src.domain_particles, dst.domain_particles);
        append_particles<ParticleType::LevelGhost>(src.level_ghost_particles,
                                                   dst.level_ghost_particles);
        return ions;
    }
}


// GPU: the ThreadedStreamLauncher's per patch, per host function times (ns), MultiBoris::move
//  registers "move" (push, or COPY's launch) then "sync" (REF: stream wait + sync_ref)
//  a time is the interval from the function's start on its patch to its end, waits on the
//  patch's stream included: patches overlap, so these are per patch, not summed
//  per iteration: accumulated in launcher_sums, the last iteration dumped to
//  $PHARE_ASYNC_TIMES/<benchmark name>.txt (lines: patch fn ns, as dump_times)
//  after the loop (no launcher): counters <fn>_ns, the mean per patch per iteration
inline std::vector<double> launcher_sums;
inline std::size_t launcher_patches = 0;

void launcher_times(benchmark::State& state, auto const&... streamer)
{
    static std::array<std::string, 2> const names{"move", "sync"};
    if constexpr (sizeof...(streamer) == 1)
    {
        auto const& launcher = (streamer, ...);
        auto const n_fns     = launcher.fns.size();
        launcher_sums.resize(n_fns, 0);
        launcher_patches = launcher.times.size() / n_fns;
        for (std::size_t i = 0; i < launcher.times.size(); ++i)
            launcher_sums[i % n_fns] += launcher.times[i].time();

        std::string name = state.name();
        std::replace_if(name.begin(), name.end(), [](char c) { return !std::isalnum(c); }, '_');
        launcher.dump_times(detail::timings_dir_str + "/" + name + ".txt");
    }
    else
    {
        for (std::size_t f = 0; f < launcher_sums.size(); ++f)
            state.counters[(f < names.size() ? names[f] : "fn" + std::to_string(f)) + "_ns"]
                = launcher_sums[f] / (state.iterations() * launcher_patches);
        launcher_sums.clear();
    }
}


template<std::size_t interp, MultiBorisMode mode>
void old_boris(benchmark::State& state, auto& patches, std::size_t const n_particles);

template<bool is_gpu, MultiBorisMode mode, typename Interpolator_t>
void multi_boris_move(benchmark::State& state, auto& patches, auto& boxings,
                      std::size_t const n_particles);


template<std::size_t dim, auto layout_mode, auto alloc_mode, MultiBorisMode mode>
void multi_boris(benchmark::State& state)
{
    auto constexpr static opts = SimOpts{
        .dimension = dim, .interp_order = 1, .layout_mode = layout_mode, .alloc_mode = alloc_mode};
    auto constexpr static is_gpu  = alloc_mode == AllocatorMode::GPU_UNIFIED;
    auto constexpr static interp  = opts.interp_order;
    using Hybrid_t                = PHARE_Types<opts>::Hybrid;
    using GridLayout_t            = Hybrid_t::GridLayout_t;
    using ParticleArray_t         = Hybrid_t::ParticleArray_t;
    auto constexpr static em_opts = TensorFieldOptions<Hybrid_t>{};
    using Ions_t                  = UsableIons_t<opts>;
    using EM_t                    = UsableElectromag<em_opts>;
    using Patch_t                 = test::HybridPatch<Ions_t, EM_t>;
    using Interpolator_t          = Interpolator<dim, interp, /*atomic=*/is_gpu>;
    using Boxing_t                = UpdaterSelectionBoxing<GridLayout_t>;

    auto const scen = load_scenario<dim>();

    std::vector<aggregate_adapter<Patch_t>> patches;
    patches.reserve(scen.patches.size());
    std::unordered_map<std::string, Boxing_t> boxings;
    std::size_t n_particles = 0;
    for (std::size_t i = 0; i < scen.patches.size(); ++i)
    {
        auto const& [box, ppc] = scen.patches[i];
        auto const layout      = make_layout<opts>(scen.dl, box);
        auto& patch            = patches.emplace_back(layout, make_ions<opts>(layout, ppc, i));
        patch.id               = std::to_string(i);
        boxings.try_emplace(patch.id, Boxing_t{layout, {layout.AMRBox()}});
        n_particles += box.size() * ppc;
    }

    state.counters["patches"]   = patches.size();
    state.counters["particles"] = n_particles;
    state.counters["pools"]     = ThreadPool::n_pools;
    state.counters["threads"]   = ThreadPool::threads_per_pool;
    // seconds per particle per move (shown with SI prefix, e.g. 55.8n = 55.8 ns)
    state.counters["t_per_particle"] = benchmark::Counter(
        n_particles, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);

    if constexpr (layout_mode == LayoutMode::AoSMapped)
        old_boris<interp, mode>(state, patches, n_particles);
    else
        multi_boris_move<is_gpu, mode, Interpolator_t>(state, patches, boxings, n_particles);
}


// AoSMapped: not a MultiBoris layout, the serial IonUpdater's push and deposit (as master), per
//  patch: BorisPusher::move with no-op selectors (no partitioning), the interpolator to deposit
//  COPY: pushed from the patch's particles into scratch copies (made once), then deposited
//  REF:  pushed in place, not deposited
template<std::size_t interp, MultiBorisMode mode>
void old_boris_patch(auto& patch, auto& scratch)
{
    using GridLayout_t  = std::decay_t<decltype(patch.layout)>;
    auto constexpr dim  = GridLayout_t::dimension;
    auto constexpr noop = [](auto& range) { return range; };

    auto& [layout, ions, em, electromag] = patch.model.state;
    BorisPusher<dim> pusher{layout.meshSize(), dt};
    Interpolator<dim, interp> interpolator;

    std::size_t s = 0;
    for (auto& pop : ions)
        for (auto* particles : {&pop.domainParticles(), &pop.levelGhostParticles()})
        {
            auto in = makeIndexRange(*particles);
            if constexpr (mode == MultiBorisMode::COPY)
            {
                auto out = makeIndexRange(scratch[s++]);
                out      = pusher.move(in, out, electromag, pop.charge(), pop.mass(), interpolator,
                                       layout, noop, noop);
                interpolator(out, pop.particleDensity(), pop.flux(), layout);
            }
            else
                pusher.move(in, in, electromag, pop.charge(), pop.mass(), interpolator, layout,
                            noop, noop);
        }
}

template<std::size_t interp, MultiBorisMode mode>
void old_boris(benchmark::State& state, auto& patches, std::size_t const n_particles)
{
    // COPY's scratch copies, per patch: domain and level ghost particles of each population
    using ParticleArray_t
        = std::decay_t<decltype(patches[0].model.state.ions[0].domainParticles())>;
    std::vector<std::vector<ParticleArray_t>> scratch(patches.size());
    for (std::size_t i = 0; i < patches.size(); ++i)
        for (auto& pop : patches[i].model.state.ions)
        {
            scratch[i].push_back(pop.domainParticles());
            scratch[i].push_back(pop.levelGhostParticles());
        }

    for (auto _ : state)
        for (std::size_t i = 0; i < patches.size(); ++i)
            old_boris_patch<interp, mode>(patches[i], scratch[i]);

    state.SetItemsProcessed(state.iterations() * n_particles);
}


template<bool is_gpu, MultiBorisMode mode, typename Interpolator_t>
void multi_boris_move(benchmark::State& state, auto& patches, auto& boxings,
                      std::size_t const n_particles)
{
    auto const quantities = [&](int i) {
        return std::forward_as_tuple(patches[i].model.state.ions,
                                     patches[i].model.state.electromag);
    };
    auto accessor    = test::make_model_level_accessor(patches, quantities);
    using Accessor_t = decltype(accessor);

    for (auto _ : state)
    {
        if constexpr (mode == MultiBorisMode::REF)
        { // REF leaves arrays for post_move_sync (as the updater), untimed
            state.PauseTiming();
            for (auto& patch : patches)
                for (auto& pop : patch.model.state.ions)
                    post_move_sync(pop.domainParticles(), pop.patchGhostParticles(),
                                   pop.levelGhostParticles(), boxings.at(patch.id));
            state.ResumeTiming();
        }

        if constexpr (is_gpu)
        {
#if PHARE_HAVE_MKN_GPU
            mkn_xyz::MultiBoris<Accessor_t, Interpolator_t> in{dt, accessor};
            in.template move<mode>(boxings);
            in.streamer.join();
            launcher_times(state, in.streamer);
#endif
        }
        else
            core::MultiBoris<Accessor_t, Interpolator_t>{dt, accessor}.template move<mode>(boxings);
    }

    state.SetItemsProcessed(state.iterations() * n_particles); // items/s = particles/s
    if constexpr (is_gpu)
        launcher_times(state);
}

} // namespace PHARE::core::bench


using PHARE::AllocatorMode;
using PHARE::core::LayoutMode;
using PHARE::core::MultiBorisMode;
using PHARE::core::bench::multi_boris;

#ifndef PHARE_BENCH_DIM
#define PHARE_BENCH_DIM 3
#endif
auto constexpr static dim = PHARE_BENCH_DIM;

// clang-format off
#define PHARE_BENCH_LAYOUT(layout, alloc)                                                          \
    BENCHMARK_TEMPLATE(multi_boris, dim, LayoutMode::layout, AllocatorMode::alloc,                 \
                       MultiBorisMode::COPY)->Unit(benchmark::kMillisecond)->UseRealTime();        \
    BENCHMARK_TEMPLATE(multi_boris, dim, LayoutMode::layout, AllocatorMode::alloc,                 \
                       MultiBorisMode::REF)->Unit(benchmark::kMillisecond)->UseRealTime();

PHARE_BENCH_LAYOUT(AoSMapped, CPU) // not MultiBoris: the serial updater's push/deposit (master)
PHARE_BENCH_LAYOUT(AoS,       CPU)
PHARE_BENCH_LAYOUT(AoSCMTS,   CPU)
PHARE_BENCH_LAYOUT(AoSPCTS,   CPU)
PHARE_BENCH_LAYOUT(SoAVX,     CPU)
#if PHARE_HAVE_GPU
PHARE_BENCH_LAYOUT(AoSTS,     GPU_UNIFIED)
PHARE_BENCH_LAYOUT(AoSPCTS,   GPU_UNIFIED)
PHARE_BENCH_LAYOUT(SoAVX,     GPU_UNIFIED)
#endif
// clang-format on

BENCHMARK_MAIN();
