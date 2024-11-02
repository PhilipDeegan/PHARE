//  tests/core/numerics/ion_updater/test_multi_updater.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

// USE HIP_VISIBLE_DEVICES OR CUDA_VISIBLE_DEVICES env vars

#include "core/numerics/ion_updater/ion_updater_def.hpp"

// #define PHARE_UNDEF_ASSERT
#define PHARE_SKIP_MPI_IN_CORE

#include <cstddef>

#include "core/data/mkn.gpu.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/ion_updater/ion_updater_per_particle.hpp"
#include "core/numerics/ion_updater/ion_updater_multi_pc.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"


namespace PHARE::core
{


template<std::size_t dim, typename internals> // used by gtest
void PrintTo(ParticleArray<dim, internals> const& arr, std::ostream* os)
{
    // assert(arr.size());
    *os << arr;
}

auto static const bytes     = get_env_as("PHARE_GPU_BYTES", std::uint64_t{500000000});
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{5});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{5});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1039});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const shufle    = get_env_as("PHARE_UNSORTED", std::size_t{0});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
auto static const device_id = get_env_as("PHARE_GPU_DEVICE", std::size_t{0});

bool static const premain = []() {
    PHARE_WITH_MKN_GPU({
        mkn::gpu::setDevice(device_id);
        mkn::gpu::setLimitMallocHeapSize(bytes);
    })
    PHARE_WITH_PHLOP(                                     //
        PHARE_LOG_LINE_STR("cells      : " << cells);     //
        PHARE_LOG_LINE_STR("n_patches  : " << n_patches); //
        PHARE_LOG_LINE_STR("ppc        : " << ppc);       //
        PHARE_LOG_LINE_STR("seed       : " << seed);

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

template<typename IonUpdater_t>
auto construct_()
{
    PHARE::initializer::PHAREDict dict;
    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    return IonUpdater_t{dict["simulation"]["algo"]["ion_updater"]};
}

template<typename Ions, typename EM, typename GridLayout_t>
auto get_updater_for(Ions& /*ions*/, EM const& /*em*/, GridLayout_t const& /*layout*/)
{
    using Particles = typename Ions::particle_array_type;
    if constexpr (Particles::is_mapped)
        return construct_<IonUpdater<Ions, EM, GridLayout_t>>();
    else if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSPC, LayoutMode::SoAPC))
        return construct_<IonUpdaterMultiPC<Ions, EM, GridLayout_t>>();
    else
        return construct_<IonUpdaterPP<Ions, EM, GridLayout_t>>();
}


namespace detail::strings
{
    constexpr static std::string_view update = "update,";
    constexpr static std::string_view cma    = ",";
} // namespace detail::strings

template<typename Ions, typename EM, typename GridLayout_t>
void update(UpdaterMode mode, Ions& ions, EM const& em, GridLayout_t const& layout)
{
    using Particles = typename Ions::particle_array_type;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);
    get_updater_for(ions, em, layout).updatePopulations(ions, em, layout, dt, mode);
}

template<typename Patches>
void update(UpdaterMode mode, Patches& patches)
{
    using Particles = Patches::value_type::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);
    get_updater_for(*patches[0].ions, *patches[0].em, patches[0].layout)
        .updatePopulations(patches, dt, mode);
}


template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;
    auto constexpr static interp     = GridLayout_t::interp_order;

    auto ions_p = std::make_shared<UsableIons_t<Particles_t, interp>>(layout, "protons");
    auto& ions  = *ions_p;

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0);

    auto disperse = [](auto& particles) {
        delta_disperse(particles.domain_particles, seed);
        delta_disperse(particles.patch_ghost_particles, seed);
        if (shufle > 0)
        {
            shuffle(particles.domain_particles, seed);
            shuffle(particles.patch_ghost_particles, seed);
        }
    };

    auto add_particles = [&](auto& particles) {
        add_particles_in(particles.domain_particles, layout.AMRBox(), ppc);
        add_ghost_particles(particles.patch_ghost_particles, layout.AMRBox(), ppc,
                            GridLayout_t::nbrParticleGhosts());
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);

    return ions_p;
}



template<typename Particles_t, typename GridLayout_t, typename Ions>
auto from_ions(GridLayout_t const& layout, Ions const& from)
{
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    // auto constexpr static dim        = GridLayout_t::dimension;
    auto constexpr static interp = GridLayout_t::interp_order;
    auto ions_p = std::make_shared<UsableIons_t<Particles_t, interp>>(layout, "protons");
    auto& ions  = *ions_p;
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0);

    auto _add_particles_from = [&]<auto type>(auto& src, auto& dst) {
        ParticleArrayService::reserve_ppc_in<type>(dst, ppc);
        append<type>(src, dst, layout);
    };

    _add_particles_from.template operator()<ParticleType::Domain>(
        from.populations[0].particles.domain_particles,
        ions.populations[0].particles.domain_particles);
    _add_particles_from.template operator()<ParticleType::Ghost>(
        from.populations[0].particles.patch_ghost_particles,
        ions.populations[0].particles.patch_ghost_particles);
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);
    return ions_p;
}

template<auto updater_mode, typename Patches>
auto& evolve(Patches& patches)
{
    update(updater_mode, patches);
    return patches;
}

template<auto updater_mode, typename Ions, typename GridLayout_t>
auto& evolve(Ions& ions, GridLayout_t const& layout)
{
    using Particles_t                = typename Ions::particle_array_type;
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;

    UsableElectromag<dim, alloc_mode> em{layout};
    update(updater_mode, *ions, *em, layout);
    assert(ions.populations[0].particles.domain_particles.size() > 0);
    return ions;
}


template<typename GridLayout_t, typename P0, typename P1>
void compare_particles(GridLayout_t const& layout, P0& ref, P1& cmp_, std::size_t ghosts = 0)
{
    auto const box = grow(layout.AMRBox(), ghosts);

    sort(cmp_, box);
    auto const cmp = ParticleArrayService::convert_to<LayoutMode::AoS>(cmp_, layout);
    sort(ref, box);

    EXPECT_EQ(ref.size(), cmp.size());

    auto report = ref == cmp;
    if (!report)
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays: " << P0::id() << " vs " << P1::id());
        PHARE_LOG_LINE_STR("results: " << report.why());
    }

    EXPECT_EQ(ref, cmp);
}

template<bool GHOSTS = false, typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    // PHARE_LOG_LINE_STR(ref.populations[0].particles.domain_particles.begin().copy());
    // PHARE_LOG_LINE_STR(cmp.populations[0].particles.domain_particles.begin().copy());

    compare_particles(layout, ref.populations[0].particles.domain_particles,
                      cmp.populations[0].particles.domain_particles);

    if constexpr (GHOSTS)
        compare_particles(layout, ref.populations[0].particles.patch_ghost_particles,
                          cmp.populations[0].particles.patch_ghost_particles);

    EXPECT_TRUE(ref.populations[0].F.isclose(cmp.populations[0].F));     //
    EXPECT_TRUE(ref.populations[0].rho.isclose(cmp.populations[0].rho)); //
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode, std::uint8_t _impl,
         auto _updater_mode>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), PHARE::AllocatorMode>);
    auto constexpr static dim          = _dim;
    auto constexpr static layout_mode  = _layout_mode;
    auto constexpr static alloc_mode   = _alloc_mode;
    auto constexpr static impl         = _impl;
    auto constexpr static updater_mode = _updater_mode;
};

template<typename Ions, std::size_t dim, std::size_t interp, auto updater_mode>
struct DefaultIons
{
    using GridLayout_t = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;

    GridLayout_t const layout{cells};
    std::shared_ptr<Ions> init = make_ions<typename Ions::particle_array_type>(layout);
    std::shared_ptr<Ions> ions = make_ions<typename Ions::particle_array_type>(layout);

    DefaultIons() { evolve<updater_mode>(*ions, *layout); }

    static auto& I()
    {
        static DefaultIons self;
        return self;
    }
};

template<typename Param>
struct MultiPatchIonUpdaterTest : public ::testing::Test
{
    auto constexpr static dim          = Param::dim;
    auto constexpr static interp       = 1;
    auto constexpr static layout_mode  = Param::layout_mode;
    auto constexpr static alloc_mode   = Param::alloc_mode;
    auto constexpr static impl         = Param::impl;
    auto constexpr static updater_mode = Param::updater_mode;


    using GridLayout_t       = TestGridLayout<typename PHARE_Types<dim, interp>::GridLayout_t>;
    using RefParticleArray_t = AoSMappedParticleArray<dim>;
    using CmpParticleArray_t = ParticleArray<
        dim, ParticleArrayInternals<dim, layout_mode, StorageMode::VECTOR, alloc_mode, impl>>;

    using RefIons_t = UsableIons_t<RefParticleArray_t, interp>;
    using CmpIons_t = UsableIons_t<CmpParticleArray_t, interp>;
    using DefIons   = DefaultIons<RefIons_t, dim, interp, Param::updater_mode>;

    GridLayout_t const& layout           = DefIons::I().layout;
    std::shared_ptr<RefIons_t>& ref_ions = DefIons::I().ions;


    MultiPatchIonUpdaterTest()
    {
        patches.resize(n_patches);
        // uncomment to do init value check

        // auto ref = make_ions<RefParticleArray_t>(layout);
        // auto cmp = from_ions<CmpParticleArray_t>(layout, *ref);
        // compare(*this->layout, *ref, *cmp); // pre update check
    }

    struct Patch
    {
        using ParticleArray_t    = CmpParticleArray_t;
        using UsableElectromag_t = UsableElectromag<dim, alloc_mode>;
        using Electromag_t       = UsableElectromag_t::Super;
        using GridLayout_t       = MultiPatchIonUpdaterTest<Param>::GridLayout_t;

        GridLayout_t layout             = DefIons::I().layout;
        std::shared_ptr<CmpIons_t> ions = from_ions<CmpParticleArray_t>(layout, *DefIons::I().init);
        UsableElectromag_t em{layout};

        Patch() { mkn::gpu::print_gpu_mem_used(); }
    };

    std::vector<Patch> patches{};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

PHARE_WITH_MKN_GPU(

    // TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::domain_only>,
    // TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::domain_only>,

    TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>,
    TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>

)

>;
// clang-format on

TYPED_TEST_SUITE(MultiPatchIonUpdaterTest, Permutations_t, );

template<typename MultiPatchIonUpdaterTest_t>
auto run(MultiPatchIonUpdaterTest_t& self)
{
    evolve<MultiPatchIonUpdaterTest_t::updater_mode>(self.patches);
    if (do_cmp)
        for (auto const& patch : self.patches)
            compare(*self.layout, *self.ref_ions, *patch.ions);
}


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
