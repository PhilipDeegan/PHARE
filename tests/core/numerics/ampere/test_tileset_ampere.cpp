//  tests/core/numerics/ampere/test_tileset_ampere.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

// USE HIP_VISIBLE_DEVICES OR CUDA_VISIBLE_DEVICES env vars
// #define PHARE_UNDEF_ASSERT

#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/numerics/ampere/ampere_transformer.hpp"
#include "core/data/field/field_tiles_gpu.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>
#include <cstddef>
#include <algorithm>

namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{4});

bool static const premain = []() {
    ThreadPool::threads_per_pool = 4;
    return true;
}();

auto& pool = *ThreadPool::INSTANCE().thread_pools[0];


template<typename Patch>
void ref_do(Patch& patch)
{
    AmpereTransformer<>{}(patch.layout, patch.em.B, patch.J);
}

template<typename Patch>
void cmp_do(Patch& patch)
{
    AmpereTransformer<> const ampere{};

    if constexpr (Patch::alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        ampere(patch.layout, patch.em.B, patch.J); // one launch, block per tile
        return;
    }

    // a pool task per tile, as the level transformer
    for (std::size_t ti = 0; ti < tile_count(patch.J); ++ti)
        pool.detach_task([&, ti]() { ampere.on_tile(ti, patch.em.B, patch.J); });
    pool.wait();
}



template<auto alloc_mode, typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using enum LayoutMode;
    using enum AllocatorMode;
    auto constexpr static n_components = 3;

    double diff = 1e-15;

    // if constexpr (alloc_mode == GPU_UNIFIED)
    //     diff *= 1e3; // atomic no order guaranteed

    auto const& patch_box = layout.AMRBox();
    for (std::size_t c = 0; c < n_components; ++c)
    {
        {
            auto const eq = compare_fields(ref.em.B[c], cmp.em.B[c]);
            EXPECT_TRUE(eq) << "Failure for B New: " << eq.why();
        }
        {
            auto const eq = compare_fields(ref.em.E[c], cmp.em.E[c]);
            EXPECT_TRUE(eq) << "Failure for E New: " << eq.why();
        }
        auto const eq = compare_fields(ref.J[c], cmp.J[c]);
        EXPECT_TRUE(eq) << "Failure for J: " << eq.why();
    }
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};


template<typename Param>
struct AmpereTileTest : public ::testing::Test
{
    struct TileSetPatch
    {
        auto constexpr static dim         = Param::dim;
        auto constexpr static interp      = 1;
        auto constexpr static layout_mode = Param::layout_mode;
        auto constexpr static alloc_mode  = Param::alloc_mode;
        auto constexpr static opts        = SimOpts{.dimension    = dim,
                                                    .interp_order = interp,
                                                    .layout_mode  = layout_mode,
                                                    .alloc_mode   = alloc_mode};

        using PhareTypes                 = PHARE_Types<opts>;
        using GridLayout_t               = PhareTypes::Hybrid::GridLayout_t;
        using Field_t                    = PhareTypes::Hybrid::Grid_t;
        using Hybrid_t                   = PhareTypes::Hybrid;
        auto constexpr static field_opts = TensorFieldOptions<Hybrid_t>{};
        using UsableElectromag_t         = UsableElectromag<field_opts>;
        using UsableVecField_t           = UsableVecField<field_opts>;

        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J, 0};
    };

    struct ContiguousPatch
    {
        // reference (untiled) is always run on CPU, compared against the tiled/possibly-GPU patch
        auto constexpr static dim    = Param::dim;
        auto constexpr static interp = 1;
        auto constexpr static opts   = SimOpts{.dimension = dim, .interp_order = interp};

        using PhareTypes                 = PHARE_Types<opts>;
        using GridLayout_t               = PhareTypes::Hybrid::GridLayout_t;
        using Hybrid_t                   = PhareTypes::Hybrid;
        auto constexpr static field_opts = TensorFieldOptions<Hybrid_t>{};
        using UsableElectromag_t         = UsableElectromag<field_opts>;
        using Electromag_t               = UsableElectromag_t::Super;
        using Grid_t                     = PhareTypes::Hybrid::Grid_t;
        using UsableVecField_t           = UsableVecField<field_opts>;

        ContiguousPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J, 0};
    };

    using RefPatch = ContiguousPatch;
    using CmpPatch = TileSetPatch;

    AmpereTileTest() {}

    void do_compare() const { compare<CmpPatch::alloc_mode>(*cmp_layout, ref_patch, cmp_patch); }

    TestGridLayout<typename RefPatch::GridLayout_t> const ref_layout{cells};
    TestGridLayout<typename CmpPatch::GridLayout_t> const cmp_layout{cells};
    RefPatch ref_patch{ref_layout};
    CmpPatch cmp_patch{cmp_layout};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
     TestParam<3, LayoutMode::AoSCMTS, AllocatorMode::CPU>
    ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>
PHARE_WITH_GPU(
    ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED>
)
>;
// clang-format on

TYPED_TEST_SUITE(AmpereTileTest, Permutations_t, );

template<typename AmpereTileTest_t>
auto run(AmpereTileTest_t& self)
{
    ref_do(self.ref_patch);
    cmp_do(self.cmp_patch);
    self.do_compare();
}

TYPED_TEST(AmpereTileTest, dispatch)
{
    run(*this);
}

// tile buffers of tiled vecfields, every element (ghosts included)
void copy_tiles(auto& dst, auto const& src)
{
    for (std::size_t c = 0; c < dst.size(); ++c)
        for (std::size_t ti = 0; ti < dst[c]().size(); ++ti)
        {
            auto& d       = dst[c]()[ti]();
            auto const& s = src[c]()[ti]();
            std::copy(s.data(), s.data() + s.size(), d.data());
        }
}

void expect_same_tiles(auto const& a, auto const& b, std::string const& what)
{
    std::size_t n_diff = 0, n = 0;
    for (std::size_t c = 0; c < a.size(); ++c)
        for (std::size_t ti = 0; ti < a[c]().size(); ++ti)
        {
            auto const& x = a[c]()[ti]();
            auto const& y = b[c]()[ti]();
            for (std::size_t i = 0; i < x.size(); ++i, ++n)
                n_diff += !(x.data()[i] == y.data()[i]
                            or (std::isnan(x.data()[i]) and std::isnan(y.data()[i])));
        }
    EXPECT_EQ(n_diff, 0u) << what << ": " << n_diff << " of " << n << " elements differ";
}

// device sync_inner_ghosts / average against their host versions, same input, exact
TYPED_TEST(AmpereTileTest, device_ops)
{
    using Patch = TestFixture::CmpPatch;
    if constexpr (Patch::alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        auto& host = this->cmp_patch;
        Patch dev{this->cmp_layout};
        // a value per (component, tile, element): ghosts differ from their neighbours' domains
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t ti = 0; ti < host.J[c]().size(); ++ti)
            {
                auto& f = host.J[c]()[ti]();
                for (std::size_t i = 0; i < f.size(); ++i)
                    f.data()[i] = c * 1e6 + ti * 1e3 + i * 1e-3;
            }
        copy_tiles(dev.J, host.J);
        copy_tiles(dev.em.E, host.em.E);
        copy_tiles(dev.em.B, host.em.B);

        sync_inner_ghosts(host.J);
        sync_inner_ghosts_on_device(dev.J);
        expect_same_tiles(host.J, dev.J, "sync_inner_ghosts");

        { // the host sync did write ghosts, else the comparison above proves nothing
            Patch untouched{this->cmp_layout};
            copy_tiles(untouched.J, dev.J);
            for (std::size_t c = 0; c < 3; ++c)
                for (std::size_t ti = 0; ti < untouched.J[c]().size(); ++ti)
                {
                    auto& f = untouched.J[c]()[ti]();
                    for (std::size_t i = 0; i < f.size(); ++i)
                        f.data()[i] = c * 1e6 + ti * 1e3 + i * 1e-3;
                }
            std::size_t n_synced = 0;
            for (std::size_t c = 0; c < 3; ++c)
                for (std::size_t ti = 0; ti < host.J[c]().size(); ++ti)
                    for (std::size_t i = 0; i < host.J[c]()[ti]().size(); ++i)
                        n_synced += host.J[c]()[ti]().data()[i]
                                    != untouched.J[c]()[ti]().data()[i];
            if (host.J[0]().size() > 1)
            {
                EXPECT_GT(n_synced, 0u);
            }
        }

        for (std::size_t c = 0; c < 3; ++c) // host reference
            for (std::size_t ti = 0; ti < host.em.B[c]().size(); ++ti)
            {
                auto const& x0 = host.em.E[c]()[ti]();
                auto const& x1 = host.J[c]()[ti]();
                auto& a        = host.em.B[c]()[ti]();
                for (std::size_t i = 0; i < a.size(); ++i)
                    a.data()[i] = (x0.data()[i] + x1.data()[i]) * .5;
            }
        core::average(dev.em.E, dev.J, dev.em.B); // GPU_UNIFIED: device jobs, see algorithm.hpp
        expect_same_tiles(host.em.B, dev.em.B, "average");
    }
}

} // namespace PHARE::core

int main(int argc, char** argv)
{
    // assert(phlop::ScopeTimerMan::INSTANCE().active);
    ::testing::InitGoogleTest(&argc, argv);
    auto r = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(if (phlop::scope_timer().active) phlop::shutdown_scope_timer();)
    return r;
}
