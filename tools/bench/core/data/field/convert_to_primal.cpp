
#include "simulator/simulator_def.hpp"

#include "tools/bench/core/bench.hpp"
#include "core/utilities/algorithm.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

// compares ways of converting E and B to all primal, as diagnostics do before writing:
//   Master   : quantity checked per index, then GridLayout::project
//   Recenter : RecenteredField (current core::convert_to_fortran_primal)
//   WholeBox : quantity checked once, then GridLayout::project over the whole box

using namespace PHARE::core;

template<std::size_t dim>
std::uint32_t constexpr cells = dim == 1   ? 100000
                                : dim == 2 ? 300
                                           : 50;

auto box_views(auto& dst, auto const& layout)
{
    bool static constexpr c_ordering = false;

    auto const lcl_box = layout.AMRToLocal(layout.AMRBoxFor(dst));
    auto lcl_zero_box  = lcl_box;
    lcl_zero_box.lower -= lcl_box.lower;
    lcl_zero_box.upper -= lcl_box.lower;
    return std::make_tuple(lcl_box, lcl_zero_box,
                           make_array_view<c_ordering>(dst.data(), *lcl_box.shape()));
}

struct Master
{
    static auto convert_to_primal(auto const& src, auto const& layout, auto const lix,
                                  HybridQuantity::Scalar const qty)
    {
        using PQ         = HybridQuantity::Scalar;
        using GridLayout = std::decay_t<decltype(layout)>;

        if (qty == PQ::Bx)
            return GridLayout::template project<GridLayout::BxToMoments>(src, lix);
        else if (qty == PQ::By)
            return GridLayout::template project<GridLayout::ByToMoments>(src, lix);
        else if (qty == PQ::Bz)
            return GridLayout::template project<GridLayout::BzToMoments>(src, lix);
        else if (qty == PQ::Ex)
            return GridLayout::template project<GridLayout::ExToMoments>(src, lix);
        else if (qty == PQ::Ey)
            return GridLayout::template project<GridLayout::EyToMoments>(src, lix);
        else if (qty == PQ::Ez)
            return GridLayout::template project<GridLayout::EzToMoments>(src, lix);

        throw std::runtime_error("Quantity not supported for conversion to primal.");
    }

    static void convert_to_fortran_primal(auto& dst, auto const& src, auto const& layout)
    {
        auto [lcl_box, lcl_zero_box, dst_view] = box_views(dst, layout);
        for_each_point_zip(lcl_box, lcl_zero_box, [&](auto const& lix, auto const& lix0) {
            dst_view(lix0) = convert_to_primal(src, layout, lix, src.physicalQuantity());
        });
    }
};

struct Recenter
{
    static void convert_to_fortran_primal(auto& dst, auto const& src, auto const& layout)
    {
        PHARE::core::convert_to_fortran_primal(dst, src, layout);
    }
};

struct WholeBox
{
    template<auto projector>
    static void project(auto& dst, auto const& src, auto const& layout)
    {
        using GridLayout = std::decay_t<decltype(layout)>;

        auto [lcl_box, lcl_zero_box, dst_view] = box_views(dst, layout);
        for_each_point_zip(lcl_box, lcl_zero_box, [&](auto const& lix, auto const& lix0) {
            dst_view(lix0) = GridLayout::template project<projector>(src, lix);
        });
    }

    static void convert_to_fortran_primal(auto& dst, auto const& src, auto const& layout)
    {
        using PQ         = HybridQuantity::Scalar;
        using GridLayout = std::decay_t<decltype(layout)>;

        auto const qty = src.physicalQuantity();
        if (qty == PQ::Bx)
            project<GridLayout::BxToMoments>(dst, src, layout);
        else if (qty == PQ::By)
            project<GridLayout::ByToMoments>(dst, src, layout);
        else if (qty == PQ::Bz)
            project<GridLayout::BzToMoments>(dst, src, layout);
        else if (qty == PQ::Ex)
            project<GridLayout::ExToMoments>(dst, src, layout);
        else if (qty == PQ::Ey)
            project<GridLayout::EyToMoments>(dst, src, layout);
        else if (qty == PQ::Ez)
            project<GridLayout::EzToMoments>(dst, src, layout);
        else
            throw std::runtime_error("Quantity not supported for conversion to primal.");
    }
};


template<std::size_t dim, typename Impl>
void convert_to_primal(benchmark::State& state)
{
    auto static constexpr opts = PHARE::SimOpts{dim, 1};
    using GridLayout_t         = TestGridLayout<typename PHARE_Types<opts>::Hybrid::GridLayout_t>;

    GridLayout_t layout{cells<dim>};
    UsableElectromag<dim> em{layout};
    UsableVecField<dim> dst{"dst", layout, HybridQuantity::Vector::V};
    for (auto* vf : {&em.E, &em.B})
        for (auto& c : *vf)
            std::fill(c.begin(), c.end(), 1.);

    while (state.KeepRunning())
        for (auto const* vf : {&em.E, &em.B})
            for (std::size_t ci = 0; ci < 3; ++ci)
            {
                Impl::convert_to_fortran_primal(dst[ci], (*vf)[ci], layout);
                benchmark::DoNotOptimize(dst[ci].data());
            }
}

BENCHMARK_TEMPLATE(convert_to_primal, 1, Master)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 1, Recenter)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 1, WholeBox)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(convert_to_primal, 2, Master)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 2, Recenter)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 2, WholeBox)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(convert_to_primal, 3, Master)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 3, Recenter)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(convert_to_primal, 3, WholeBox)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
