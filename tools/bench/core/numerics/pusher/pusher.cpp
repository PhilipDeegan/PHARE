
#include "benchmark/benchmark.h"

#include "phare_core.h"
#include "core/numerics/pusher/boris.h"
#include "core/numerics/ion_updater/ion_updater.h"

static constexpr std::uint32_t min_icell = 5; // first physical index;

template<typename Float, std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim, Float>,
                                 typename PHARE::core::HybridQuantity::Scalar>;

template<typename Float, std::size_t dim>
PHARE::core::Particle<Float, dim> particle()
{
    return {//
            /*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ PHARE::core::ConstArray<int, dim>(min_icell),
            /*.delta  = */ PHARE::core::ConstArray<float, dim>(.01),
            /*.v      = */ {{0, 10., 0}}};
}

template<typename ParticleArray>
void disperse(ParticleArray& particles, std::size_t cells)
{
    std::random_device rd;
    std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 gen{seed_seq};
    std::uniform_int_distribution<> distrib(min_icell, cells);

    for (auto& particle : particles)
        for (std::size_t i = 0; i < ParticleArray::dimension; i++)
            particle.iCell[i] = distrib(gen);
}

template<typename GridLayout, typename Quantity, typename Float = typename GridLayout::Float,
         std::size_t dim = GridLayout::dimension>
Field<Float, dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<Float, dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}

template<typename Float, std::size_t dim, std::size_t interp, bool dispersal = true>
void push(benchmark::State& state)
{
    constexpr Float one           = 1;
    constexpr std::uint32_t cells = 65;
    constexpr std::uint32_t parts = 1e8;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp, Float>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp, Float>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using Electromag_t      = typename PHARE_Types::Electromag_t;
    using GridLayout_t      = typename PHARE_Types::GridLayout_t;
    using Field_t           = typename PHARE_Types::Field_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using PartIterator      = typename ParticleArray::iterator;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, PartIterator, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;

    Interpolator interpolator;
    ParticleArray domainParticles{parts, particle<Float, dim>()};
    if constexpr (dispersal)
        disperse(domainParticles, cells);
    ParticleArray tmpDomain{domainParticles.size(), particle<Float, dim>()};

    auto rangeIn  = PHARE::core::makeRange(domainParticles);
    auto rangeOut = PHARE::core::makeRange(tmpDomain);

    auto meshSize = PHARE::core::ConstArray<Float, dim>(one / cells);
    auto nCells   = PHARE::core::ConstArray<std::uint32_t, dim>(cells);
    auto origin   = PHARE::core::Point<Float, dim>{PHARE::core::ConstArray<Float, dim>(0)};
    GridLayout_t layout{meshSize, nCells, origin};

    Field_t bx = field("Bx", PHARE::core::HybridQuantity::Scalar::Bx, layout);
    Field_t by = field("By", PHARE::core::HybridQuantity::Scalar::By, layout);
    Field_t bz = field("Bz", PHARE::core::HybridQuantity::Scalar::Bz, layout);

    Field_t ex = field("Ex", PHARE::core::HybridQuantity::Scalar::Ex, layout);
    Field_t ey = field("Ey", PHARE::core::HybridQuantity::Scalar::Ey, layout);
    Field_t ez = field("Ez", PHARE::core::HybridQuantity::Scalar::Ez, layout);

    Electromag_t emFields{std::string{"EM"}};
    emFields.B.setBuffer("EM_B_x", &bx);
    emFields.B.setBuffer("EM_B_y", &by);
    emFields.B.setBuffer("EM_B_z", &bz);
    emFields.E.setBuffer("EM_E_x", &ex);
    emFields.E.setBuffer("EM_E_y", &ey);
    emFields.E.setBuffer("EM_E_z", &ez);

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    for (auto _ : state)
    {
        pusher.move(
            /*ParticleRange const&*/ rangeIn, /*ParticleRange&*/ rangeOut,
            /*Electromag const&*/ emFields, /*double mass*/ 1, /*Interpolator&*/ interpolator,
            /*ParticleSelector const&*/ [](auto const& /*part*/) { return true; },
            /*GridLayout const&*/ layout);
    }
}
BENCHMARK_TEMPLATE(push, double, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, double, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, double, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, double, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);


BENCHMARK_TEMPLATE(push, float, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, float, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, float, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, float, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
