//  tests/core/numerics/ampere/test_tileset_ampere.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

// USE HIP_VISIBLE_DEVICES OR CUDA_VISIBLE_DEVICES env vars
// #define PHARE_UNDEF_ASSERT
#include "core/utilities/types.hpp"
#define PHARE_SKIP_MPI_IN_CORE

#include "core/def/phare_config.hpp"

#include "initializer/data_provider.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>

#include "bstp/include/BS_thread_pool.hpp"

namespace PHARE::core
{
// COMPILE TIME SWITCHES
std::size_t static constexpr tile_size = 4;

// RUNTIME ENV VAR OVERRIDES
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
static ::BS::thread_pool pool{4};

bool static const premain = []() {
    PHARE_WITH_PHLOP(                                     //
        PHARE_LOG_LINE_STR("cells      : " << cells);     //
        PHARE_LOG_LINE_STR("n_patches  : " << n_patches); //

        using namespace PHARE; //
        using namespace std::literals;
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::threaded::ScopeTimerMan::INSTANCE()
                .file_name(".phare_times.0.txt")
                // .force_strings()
                // .headers("fn"s, "dim"s, "layout"s, "alloc"s, "storage"s, "time"s)
                .init(); //
    )
    return true;
}();



namespace detail::strings
{
    constexpr static std::string_view test    = "test,";
    constexpr static std::string_view compare = "compare,";
    constexpr static std::string_view cma     = ",";
} // namespace detail::strings


template<typename Patches, typename Patch = Patches::value_type>
void ref_do(Patches& patches)
{
    using GridLayout_t = Patch::GridLayout_t;

    for (auto& patch : patches)
        Ampere_ref<GridLayout_t>{patch.layout}(patch.em.B, patch.J);
}

template<typename Patches, typename Patch = Patches::value_type>
void cmp_do(Patches& patches)
{
    using GridLayout_t = Patch::GridLayout_t::Super;
    using GridTile_t   = Patch::Field_t::value_type;
    using Field_vt     = GridTile_t::value_type::field_type;
    using VecField_vt  = basic::TensorField<Field_vt>;

    auto constexpr tile_accessor = [](auto& ts, auto const ti) { return *ts[ti](); };

    for (std::size_t pi = 0; pi < patches.size(); ++pi)
        pool.detach_task([&, idx = pi]() {
            auto& patch        = patches[idx];
            auto const n_tiles = patch.J[0].size();
            for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
            {
                auto const B = patch.em.B.template as<VecField_vt>(tile_accessor, ti);
                auto J       = patch.J.template as<VecField_vt>(tile_accessor, ti);
                pool.detach_task([=, layout = patch.em.B[0][ti].layout()]() mutable {
                    Ampere_ref<GridLayout_t>{layout}(B, J);
                });
            }
        });
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
    std::size_t non_zero  = 0;

    for (std::size_t c = 0; c < n_components; ++c)
    {
        auto const& ref_Jxyz = ref.J[c];
        auto const& cmp_Jxyz = cmp.J[c];

        for (auto const& tile : cmp_Jxyz)
        {
            auto const& tile_field = tile();
            for (auto const& tix : *tile)
            {
                auto const tlix = tix - tile.lower;
                auto const plix = tix - patch_box.lower;
                EXPECT_TRUE(float_equals(tile_field(tlix), ref_Jxyz(plix), diff));
                non_zero += tile_field(tlix) == 0 ? 0 : 1;
            }
        }
    }

    EXPECT_TRUE(non_zero > 0);
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode, std::uint8_t _impl>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
    auto constexpr static impl        = _impl;
};

auto inline static field_dict(std::string const& name)
{
    initializer::PHAREDict dict;
    dict["tile_size"] = tile_size;
    dict["name"]      = name;
    return dict;
}

template<typename GridLayout_t, typename Grid_t>
struct UsableTiledElectromag
{
    using VecField_t = basic::TensorField<Grid_t>;

    UsableTiledElectromag(GridLayout_t const& layout)
        : E{{
              Grid_t{field_dict("E_x"), *layout, HybridQuantity::Scalar::Ex},
              Grid_t{field_dict("E_y"), *layout, HybridQuantity::Scalar::Ey},
              Grid_t{field_dict("E_z"), *layout, HybridQuantity::Scalar::Ez},
          }}
        , B{{
              Grid_t{field_dict("B_x"), *layout, HybridQuantity::Scalar::Bx},
              Grid_t{field_dict("B_y"), *layout, HybridQuantity::Scalar::By},
              Grid_t{field_dict("B_z"), *layout, HybridQuantity::Scalar::Bz},
          }}
    {
    }

    VecField_t E;
    VecField_t B;
};

template<typename Param>
struct AmpereTileTest : public ::testing::Test
{
    auto constexpr static c_order     = true;
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static impl        = Param::impl;

    template<typename T, auto am>
    using nd_array_t   = NdArrayVector<dim, T, c_order, am>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;

    template<auto am>
    struct TileSetPatch
    {
        using GridLayout_t       = AmpereTileTest<Param>::GridLayout_t;
        using Field_t            = GridTileSet<typename GridLayout_t::Super, nd_array_t<double, am>,
                                               HybridQuantity::Scalar>;
        using UsableElectromag_t = UsableTiledElectromag<GridLayout_t, Field_t>;
        using VecField_t         = basic::TensorField<Field_t>;


        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            // B initialized later
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout};
        VecField_t J{{
            Field_t{field_dict("J_x"), *layout, HybridQuantity::Scalar::Jx},
            Field_t{field_dict("J_y"), *layout, HybridQuantity::Scalar::Jy},
            Field_t{field_dict("J_z"), *layout, HybridQuantity::Scalar::Jz},
        }};
    };

    template<auto am = AllocatorMode::CPU>
    struct ContiguousPatch
    {
        using GridLayout_t       = AmpereTileTest<Param>::GridLayout_t;
        using UsableElectromag_t = UsableElectromag<dim, am>;
        using Electromag_t       = UsableElectromag_t::Super;
        using Grid_t             = Grid<nd_array_t<double, am>, HybridQuantity::Scalar>;
        using UsableVecField_t   = UsableVecField<dim, am>;

        ContiguousPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            default_em_init(em, layout);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout};
        UsableVecField_t J{"J", *layout, HybridQuantity::Vector::J};
    };

    using RefPatch = ContiguousPatch<>;
    using CmpPatch = TileSetPatch<alloc_mode>;


    AmpereTileTest() {}


    void init_cmp_em(auto& cmp_patch)
    {
        auto constexpr static n_components = 3;
        auto const& ref_B                  = ref_patches[0].em.B;
        auto const& patch_box              = ref_patches[0].layout.AMRBox();

        for (std::size_t c = 0; c < n_components; ++c)
        {
            auto const& ref_Bxyz = ref_B[c];
            for (auto& tile : cmp_patch.em.B[c])
            {
                auto& tile_field = tile();
                for (auto const& tix : *tile)
                    tile_field(tix - tile.lower) = ref_Bxyz(tix - patch_box.lower);
            }
        }
    }

    void setup()
    {
        cmp_patches.reserve(n_patches);

        ref_patches.reserve(1);
        auto& ref = ref_patches.emplace_back(layout);
        ref_do(ref_patches);

        for (std::size_t i = 0; i < n_patches; i++)
            init_cmp_em(cmp_patches.emplace_back(layout));
    }


    void init_check() const
    {
        bool static constexpr init_value_check = 0;
        if constexpr (init_value_check)
        {
            // auto ref = make_ions<RefParticleArray_t>(layout);
            // auto cmp = from_ions<CmpParticleArray_t>(layout, *ref);
            // compare(*this->layout, *ref, *cmp);
        }
    }

    void do_compare() const
    {
        for (std::size_t i = 0; i < n_patches; i++)
            compare<alloc_mode>(*layout, ref_patches[0], cmp_patches[i]);
    }

    GridLayout_t const layout{cells};
    std::vector<aggregate_adapter<RefPatch>> ref_patches{};
    std::vector<aggregate_adapter<CmpPatch>> cmp_patches{};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU, 2>
    // ,TestParam<3, LayoutMode::SoATS, AllocatorMode::CPU, 2>

PHARE_WITH_MKN_GPU(
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, 2>
    // ,TestParam<3, LayoutMode::SoATS, AllocatorMode::GPU_UNIFIED, 2>
)

>;
// clang-format on

TYPED_TEST_SUITE(AmpereTileTest, Permutations_t, );

template<typename AmpereTileTest_t>
auto run(AmpereTileTest_t& self)
{
    self.setup();
    cmp_do(self.cmp_patches);
    if (do_cmp)
        self.do_compare();
}


TYPED_TEST(AmpereTileTest, dispatch)
{
    run(*this);
}




} // namespace PHARE::core


int main(int argc, char** argv)
{
    // assert(phlop::ScopeTimerMan::INSTANCE().active);
    ::testing::InitGoogleTest(&argc, argv);
    auto r = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
    return r;
}
