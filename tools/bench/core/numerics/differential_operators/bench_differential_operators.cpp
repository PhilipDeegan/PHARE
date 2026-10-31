#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "benchmark/benchmark.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace PHARE::core::bench
{

constexpr SimOpts mhdOpts(std::size_t dim)
{
    return SimOpts{dim,
                   0,
                   0,
                   MHDOpts::ReconstructionType::Constant,
                   MHDOpts::SlopeLimiterType::None,
                   MHDOpts::RiemannSolverType::Rusanov};
}

template<std::size_t dim>
struct Hybrid
{
    using Model = PHARE_Types<SimOpts{dim, 1}>::Hybrid;
    using Qty   = HybridQuantity;
};

template<std::size_t dim>
struct MHD
{
    using Model = PHARE_Types<mhdOpts(dim)>::MHD;
    using Qty   = MHDQuantity;
};

template<typename Model>
struct OwnedVecField
{
    template<typename Layout, typename Qty>
    OwnedVecField(Layout const& layout, Qty qty, double phase)
        : view{"v", qty}
        , grids{makeGrid(layout, 0), makeGrid(layout, 1), makeGrid(layout, 2)}
    {
        view.setBuffer(&grids);
        for (std::size_t c = 0; c < 3; ++c)
        {
            std::size_t i = 0;
            for (auto& v : grids[c])
                v = std::sin(0.01 * static_cast<double>(i++) + phase + static_cast<double>(c));
        }
    }

    double first() const { return *grids[0].begin(); }

    double sum() const
    {
        double s = 0;
        for (auto const& grid : grids)
            for (auto const& v : grid)
                s += v;
        return s;
    }

    Model::VecField_t view;
    std::array<typename Model::Grid_t, 3> grids;

private:
    template<typename Layout>
    auto makeGrid(Layout const& layout, std::size_t i) const
    {
        return typename Model::Grid_t{view[i].name(), layout, view[i].physicalQuantity(),
                                      std::numeric_limits<double>::quiet_NaN()};
    }
};

template<typename Layout>
Layout makeLayout(std::uint32_t const n)
{
    constexpr auto dim = Layout::dimension;
    std::array<double, dim> dl;
    std::array<std::uint32_t, dim> cells;
    dl.fill(0.1);
    cells.fill(n);
    return Layout{dl, cells, Point<double, dim>{}};
}

template<typename Layout>
std::int64_t updatedPoints(std::uint32_t const n)
{
    return static_cast<std::int64_t>(3. * std::pow(static_cast<double>(n), Layout::dimension));
}

template<typename Param>
void faraday(benchmark::State& state)
{
    using Model  = Param::Model;
    using Qty    = Param::Qty;
    using Layout = Model::GridLayout_t;

    auto const n      = static_cast<std::uint32_t>(state.range(0));
    auto const layout = makeLayout<Layout>(n);

    OwnedVecField<Model> B{layout, Qty::Vector::B, 0.1};
    OwnedVecField<Model> E{layout, Qty::Vector::E, 0.2};
    OwnedVecField<Model> Bnew{layout, Qty::Vector::B, 0.3};

    Faraday<Layout> op{layout};

    for (auto _ : state)
    {
        op(B.view, E.view, Bnew.view, 0.01);
        benchmark::DoNotOptimize(Bnew.first());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * updatedPoints<Layout>(n));
}

// Ampere as before FieldBoxValueAccessor: J indexed per point rather than span-wise
template<typename Layout>
struct AmperePerPoint
{
    template<typename VecField>
    void operator()(VecField const& B, VecField& J) const
    {
        auto const amount = ConstArray<std::uint32_t, Layout::dimension>(1);
        auto const curlB  = curl(layout, B);

        auto constexpr JEq
            = [](auto const& ijk, auto& Jc, auto const& curlBc) { Jc(ijk) = curlBc(ijk); };

        layout.evalOnShrinkedGhostBox(J(Component::X), amount, JEq, J(Component::X), curlB.x);
        layout.evalOnShrinkedGhostBox(J(Component::Y), amount, JEq, J(Component::Y), curlB.y);
        layout.evalOnShrinkedGhostBox(J(Component::Z), amount, JEq, J(Component::Z), curlB.z);
    }

    Layout const& layout;
};

enum class AmpereImpl { per_point, spans };

template<typename Param, AmpereImpl impl>
void ampere(benchmark::State& state)
{
    using Model  = Param::Model;
    using Qty    = Param::Qty;
    using Layout = Model::GridLayout_t;

    auto const n      = static_cast<std::uint32_t>(state.range(0));
    auto const layout = makeLayout<Layout>(n);

    OwnedVecField<Model> B{layout, Qty::Vector::B, 0.1};
    OwnedVecField<Model> J{layout, Qty::Vector::J, 0.2};

    AmperePerPoint<Layout> const per_point{layout};
    Ampere<Layout> spans{layout};

    auto const run = [&](auto& Jv) {
        if constexpr (impl == AmpereImpl::per_point)
            per_point(B.view, Jv);
        else
            spans(B.view, Jv);
    };

    { // regression: both implementations must produce the same J, bit for bit
        OwnedVecField<Model> Jref{layout, Qty::Vector::J, 0.2};
        per_point(B.view, Jref.view);
        run(J.view);
        if (J.sum() != Jref.sum())
        {
            state.SkipWithError("ampere: sum(J) differs from the per point implementation");
            return;
        }
        state.counters["sum_J"] = J.sum();
    }

    for (auto _ : state)
    {
        run(J.view);
        benchmark::DoNotOptimize(J.first());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * updatedPoints<Layout>(n));
}

#define PHARE_DIFF_OPS_CASES(name)                                                                 \
    BENCHMARK_TEMPLATE(name, Hybrid<3>)->Name(#name "/hybrid3d")->Arg(32)->Arg(128);               \
    BENCHMARK_TEMPLATE(name, MHD<3>)->Name(#name "/mhd3d")->Arg(32)->Arg(128);                     \
    BENCHMARK_TEMPLATE(name, Hybrid<2>)->Name(#name "/hybrid2d")->Arg(256)->Arg(1024);             \
    BENCHMARK_TEMPLATE(name, MHD<2>)->Name(#name "/mhd2d")->Arg(256)->Arg(1024)

#define PHARE_AMPERE_CASES(dim, nA, nB)                                                            \
    BENCHMARK_TEMPLATE(ampere, Hybrid<dim>, AmpereImpl::per_point)                                 \
        ->Name("ampere/hybrid" #dim "d/per_point")                                                 \
        ->Arg(nA)                                                                                  \
        ->Arg(nB);                                                                                 \
    BENCHMARK_TEMPLATE(ampere, Hybrid<dim>, AmpereImpl::spans)                                     \
        ->Name("ampere/hybrid" #dim "d/spans")                                                     \
        ->Arg(nA)                                                                                  \
        ->Arg(nB);                                                                                 \
    BENCHMARK_TEMPLATE(ampere, MHD<dim>, AmpereImpl::per_point)                                    \
        ->Name("ampere/mhd" #dim "d/per_point")                                                    \
        ->Arg(nA)                                                                                  \
        ->Arg(nB);                                                                                 \
    BENCHMARK_TEMPLATE(ampere, MHD<dim>, AmpereImpl::spans)                                        \
        ->Name("ampere/mhd" #dim "d/spans")                                                        \
        ->Arg(nA)                                                                                  \
        ->Arg(nB)

PHARE_DIFF_OPS_CASES(faraday);
PHARE_AMPERE_CASES(3, 32, 128);
PHARE_AMPERE_CASES(2, 256, 1024);

} // namespace PHARE::core::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
