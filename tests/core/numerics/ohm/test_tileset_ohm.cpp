//  tests/core/numerics/ohm/test_tileset_ohm.cpp
//
// #define PHARE_UNDEF_ASSERT

#include "core/utilities/types.hpp"
#define PHARE_SKIP_MPI_IN_CORE

#include "core/def/phare_config.hpp"

#include "core/numerics/ohm/ohm.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/grid/grid_tiles.hpp"
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
double static constexpr eta            = 1;
double static constexpr nu             = 1;

// RUNTIME ENV VAR OVERRIDES
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
static ::BS::thread_pool pool{4};

bool static const premain = []() { // placeholder if needed
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
    PHARE_LOG_LINE_SS("ref_do");
    using GridLayout_t = Patch::GridLayout_t;

    for (auto& patch : patches)
        Ohm_ref<GridLayout_t>{patch.layout, eta, nu}( //
            patch.n, patch.V, patch.P, patch.em.B, patch.J, patch.emNew.E);
}



template<typename Patches, typename Patch = Patches::value_type>
void cmp_do(Patches& patches)
{
    PHARE_LOG_LINE_SS("cmp_do");
    using GridLayout_t = Patch::GridLayout_t::Super;
    using GridTile_t   = Patch::Field_t::value_type;
    using Field_vt     = GridTile_t::value_type::field_type;
    using VecField_vt  = basic::TensorField<Field_vt>;

    auto constexpr tile_accessor = [](auto& ts, auto const ti) { return *ts[ti](); };

    for (std::size_t pi = 0; pi < patches.size(); ++pi)
        pool.detach_task([&, idx = pi]() {
            auto& patch        = patches[idx];
            auto const n_tiles = patch.em.B[0].size();
            for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
            {
                auto const n = *patch.n[ti]();
                auto const P = *patch.P[ti]();
                auto const V = patch.V.template as<VecField_vt>(tile_accessor, ti);
                auto const J = patch.J.template as<VecField_vt>(tile_accessor, ti);
                auto const B = patch.em.B.template as<VecField_vt>(tile_accessor, ti);
                auto ENew    = patch.emNew.E.template as<VecField_vt>(tile_accessor, ti);

                pool.detach_task([=, layout = patch.em.E[0][ti].layout()]() mutable {
                    Ohm_ref<GridLayout_t>{layout, eta, nu}(n, V, P, B, J, ENew);
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

    for (std::size_t c = 0; c < n_components; ++c)
    {
        PHARE_LOG_LINE_SS(c);
        auto const& ref_Exyz = ref.emNew.E[c];
        auto const& cmp_Exyz = cmp.emNew.E[c];

        for (auto const& tile : cmp_Exyz)
        {
            auto const& tile_field = tile();
            for (auto const& tix : *tile)
            {
                auto const tlix = tix - tile.lower;
                auto const plix = tix - patch_box.lower;
                EXPECT_TRUE(float_equals(tile_field(tlix), ref_Exyz(plix), diff));
            }
        }
    }
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

    VecField_t E, B;
};

template<typename Param>
struct OhmTileTest : public ::testing::Test
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
        using GridLayout_t       = OhmTileTest<Param>::GridLayout_t;
        using Field_t            = GridTileSet<typename GridLayout_t::Super, nd_array_t<double, am>,
                                               HybridQuantity::Scalar>;
        using UsableElectromag_t = UsableTiledElectromag<GridLayout_t, Field_t>;
        using VecField_t         = basic::TensorField<Field_t>;


        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            // EM initialized later
            J.fill(.001);
            V.fill(.001);
            n.fill(.001);
            P.fill(.001);
            em.B.fill(.001);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout}, emNew{layout};
        VecField_t J{{
            Field_t{field_dict("J_x"), *layout, HybridQuantity::Scalar::Jx},
            Field_t{field_dict("J_y"), *layout, HybridQuantity::Scalar::Jy},
            Field_t{field_dict("J_z"), *layout, HybridQuantity::Scalar::Jz},
        }};
        VecField_t V{{
            Field_t{field_dict("V_x"), *layout, HybridQuantity::Scalar::Vx},
            Field_t{field_dict("V_y"), *layout, HybridQuantity::Scalar::Vy},
            Field_t{field_dict("V_z"), *layout, HybridQuantity::Scalar::Vz},
        }};
        Field_t n{field_dict("rho"), *layout, HybridQuantity::Scalar::rho};
        Field_t P{field_dict("P"), *layout, HybridQuantity::Scalar::P};
    };

    template<auto am = AllocatorMode::CPU>
    struct ContiguousPatch
    {
        using GridLayout_t       = OhmTileTest<Param>::GridLayout_t;
        using UsableElectromag_t = UsableElectromag<dim, am>;
        using Electromag_t       = UsableElectromag_t::Super;
        using Grid_t             = Grid<nd_array_t<double, am>, HybridQuantity::Scalar>;
        using UsableVecField_t   = UsableVecField<dim, am>;

        ContiguousPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            // default_em_init(em, layout);
            J.fill(.001);
            V.fill(.001);
            n.fill(.001);
            P.fill(.001);
            em.B.fill(.001);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout, 0}, emNew{layout, 0};
        UsableVecField_t J{"J", *layout, HybridQuantity::Vector::J};
        UsableVecField_t V{"V", *layout, HybridQuantity::Vector::V};
        Grid_t n{"rho", *layout, HybridQuantity::Scalar::rho};
        Grid_t P{"P", *layout, HybridQuantity::Scalar::P};
    };



    OhmTileTest() {}


    void init_cmp_em(auto& cmp_patch)
    {
        // auto constexpr static n_components = 3;
        // auto const& patch_box              = ref_patches[0].layout.AMRBox();
        // auto const set                     = [&](auto const& ref_xyz, auto& cmp_xyz) {
        //     for (auto& tile : cmp_xyz)
        //     {
        //         auto& tile_field           = tile();
        //         auto const& tile_ghost_box = tile.layout().ghostBoxFor(tile_field);
        //         auto const& ptch_ghost_box = ref_patches[0].layout.ghostBoxFor(ref_xyz);
        //         for (auto const& tix : tile_ghost_box)
        //             tile_field(tix - tile_ghost_box.lower) = ref_xyz(tix - ptch_ghost_box.lower);
        //     }
        // };

        // for (std::size_t c = 0; c < n_components; ++c)
        // {
        //     set(ref_patches[0].em.B[c], cmp_patch.em.B[c]);
        //     set(ref_patches[0].em.E[c], cmp_patch.em.E[c]);
        // }
    }


    void setup()
    {
        cmp_patches.reserve(n_patches);

        ref_patches.reserve(1);
        ref_patches.emplace_back(layout);
        ref_do(ref_patches);

        for (std::size_t i = 0; i < n_patches; i++)
            init_cmp_em(cmp_patches.emplace_back(layout));
    }



    void do_compare() const
    {
        for (std::size_t i = 0; i < n_patches; i++)
            compare<alloc_mode>(*layout, ref_patches[0], cmp_patches[i]);
    }

    GridLayout_t const layout{cells};

    std::vector<ContiguousPatch<>> ref_patches{};
    std::vector<TileSetPatch<alloc_mode>> cmp_patches{};
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

TYPED_TEST_SUITE(OhmTileTest, Permutations_t, );

template<typename OhmTileTest_t>
auto run(OhmTileTest_t& self)
{
    self.setup();
    cmp_do(self.cmp_patches);
    if (do_cmp)
        self.do_compare();
}


TYPED_TEST(OhmTileTest, dispatch)
{
    run(*this);
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
