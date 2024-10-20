
// #include "core/data/mkn.gpu.hpp"

#include "phare_core.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_service.hpp"
#include "core/data/particles/particle_array_selector.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "gtest/gtest.h"

#define PRINT(x) std::cout << __LINE__ << " " << x << std::endl;

namespace PHARE::core
{

auto static const bytes = get_env_as("PHARE_GPU_BYTES", std::uint64_t{8000000000});
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{10});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{10});
auto static const seed  = get_env_as("PHARE_SEED", std::size_t{1039});
// auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{3});
auto static const dt = get_env_as("PHARE_TIMESTEP", double{.001});

bool static const premain = []() {
    PHARE_WITH_MKN_GPU({
        PHARE_LOG_LINE_STR("bytes: " << bytes); //
        mkn::gpu::setDevice(0);
        mkn::gpu::prinfo();
        mkn::gpu::setLimitMallocHeapSize(bytes);

        mkn::gpu::print_gpu_mem_used();
    })
    PHARE_WITH_PHLOP(                           //
        PHARE_LOG_LINE_STR("cells: " << cells); //
        PHARE_LOG_LINE_STR("ppc  : " << ppc);   //
        PHARE_LOG_LINE_STR("seed : " << seed);

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

template<std::size_t _dim, auto lm, auto am, std::uint8_t _impl = 0>
struct TestParam
{
    static_assert(std::is_same_v<decltype(lm), LayoutMode>);
    static_assert(std::is_same_v<decltype(am), PHARE::AllocatorMode>);
    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
    auto constexpr static impl        = _impl;
};


template<typename Param>
struct MultiPatchIonUpdaterTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static impl        = Param::impl;

    using GridLayout_t       = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;
    using RefParticleArray_t = AoSMappedParticleArray<dim>;
    using CmpParticleArray_t = ParticleArray<
        dim, ParticleArrayInternals<dim, layout_mode, StorageMode::VECTOR, alloc_mode, impl>>;



    MultiPatchIonUpdaterTest()
    {
        // cube of 3*3*3 cubes
        patches.reserve(3 * 3 * 3);
        int const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    int const cellx = i * cells;
                    int const celly = j * cells;
                    int const cellz = k * cells;
                    patches.emplace_back(Box<int, dim>{
                        Point{cellx, celly, cellz}, Point{cellx + off, celly + off, cellz + off}});
                }

        assert(!any_overlaps(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));
    }

    struct Patch
    {
        GridLayout_t layout;
        CmpParticleArray_t domain = make_particles<CmpParticleArray_t>(*layout);
        CmpParticleArray_t ghost  = make_particles<CmpParticleArray_t>(*layout);
    };

    bool at_periodic_border(Box<int, dim> const& box) const
    {
        auto const ghostbox = grow(box, 1);
        return domainBox * ghostbox != ghostbox;
    }

    GridLayout_t layout{cells * 3};
    Box<int, dim> const domainBox = layout.AMRBox();
    std::vector<Patch> patches{};
};



template<typename MultiPatchIonUpdaterTest_t>
auto run(MultiPatchIonUpdaterTest_t& self)
{
    auto const select = [](auto const& src, auto& dst /*, auto const& box*/) {
        auto const dst_box = grow(dst.layout.AMRBox(), 1);
        if (auto const overlap = dst_box * src.layout.AMRBox())
            ParticleArraySelector{src.domain}(*overlap, dst.ghost);
    };
    auto const select_transformed
        = [](auto const& src, auto& dst, auto const& box, auto&& transformer) {
              auto const dst_box = grow(box, 1);
              if (auto const overlap = dst_box * src.layout.AMRBox())
                  ParticleArraySelector{src.domain}(*overlap, dst.ghost, transformer);
          };

    for (auto& patch : self.patches)
    {
        add_particles_in(patch.domain, patch.layout.AMRBox(), ppc);
        PHARE_LOG_LINE_SS(patch.domain.size());
    }

    for (auto& patch : self.patches)
    {
        auto const domain_border_patch = self.at_periodic_border(patch.layout.AMRBox());

        for (auto const& src : self.patches)
        {
            if (&patch == &src)
                continue;
            select(src, patch /*, patch.layout.AMRBox()*/);
        }

        if (domain_border_patch)
        {
            // !triple transform!
            for (auto const& src : self.patches)
            {
                if (&patch == &src || !self.at_periodic_border(src.layout.AMRBox()))
                    continue;
                for_N<3>([&](auto i) {
                    int const span     = cells * 3;
                    int const shift    = patch.layout.AMRBox().lower[i] < cells ? 1 : -1;
                    auto const dst_box = shift_idx<i>(patch.layout.AMRBox(), shift * span);
                    select_transformed(src, patch, dst_box, [&](auto const& particle) {
                        auto copy    = particle;
                        copy.iCell() = (Point{copy.iCell()} + (shift * span)).toArray();
                        return copy;
                    });
                });
            }
        }
    }

    for (auto& patch : self.patches)
    {
        ParticleArrayService::template sync<2, ParticleType::Ghost>(patch.ghost);
        PHARE_LOG_LINE_STR(patch.ghost.size());
        EXPECT_EQ(9 * ppc, patch.ghost.size());
    }
}

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

PHARE_WITH_MKN_GPU(

    // TestParam<3, LayoutMode::AoSPC, AllocatorMode::CPU>,
    // TestParam<3, LayoutMode::SoAPC, AllocatorMode::CPU>,
    // TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2>,
    TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2>/*,

    TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2,UpdaterMode::all>,
    TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2,UpdaterMode::all>*/

)

>;
// clang-format on

TYPED_TEST_SUITE(MultiPatchIonUpdaterTest, Permutations_t, );

TYPED_TEST(MultiPatchIonUpdaterTest, updater_domain_only)
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
