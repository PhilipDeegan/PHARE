//  tests/core/numerics/ion_updater/test_multi_updater.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

// USE HIP_VISIBLE_DEVICES OR CUDA_VISIBLE_DEVICES env vars

// #define PHARE_UNDEF_ASSERT
#define PHARE_SKIP_MPI_IN_CORE


#include "core/utilities/kernels.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/numerics/ion_updater/ion_updater_multi_pc.hpp"
#include "core/numerics/ion_updater/ion_updater_multi_ts.hpp"
#include "core/numerics/ion_updater/ion_updater_per_particle.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"
#include <cstddef>


namespace PHARE::core
{
// COMPILE TIME SWITCHES
bool constexpr INCLUDE_PATCH_GHOST  = 0;
bool constexpr COMPARE_TO_ONE_PATCH = 1;


// RUNTIME ENV VAR OVERRIDES
auto static const bytes     = get_env_as("PHARE_GPU_BYTES", std::uint64_t{50000000});
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{4});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1009});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const shufle    = get_env_as("PHARE_UNSORTED", std::size_t{0});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
auto static const device_id = get_env_as("PHARE_GPU_DEVICE", std::size_t{0});

auto static seeder = seed; // increments over ref patches if COMPARE_TO_ONE_PATCH == 0

bool static const premain = []() {
    PHARE_WITH_MKN_GPU({
        // mkn::gpu::setDevice(device_id);
        // ::mkn::gpu::setLimitMallocHeapSize(bytes);
    })
    PHARE_WITH_PHLOP({
        PHARE_LOG_LINE_STR("cells      : " << cells);
        PHARE_LOG_LINE_STR("ppc        : " << ppc);
        PHARE_LOG_LINE_STR("particles  : " << std::pow(cells, 3) * ppc);
        PHARE_LOG_LINE_STR("n_patches  : " << n_patches);
        PHARE_LOG_LINE_STR("seed       : " << seed);
        PHARE_LOG_LINE_SS("particle MB ≈ " << std::pow(cells, 3) * ppc * 76 / 1e6);

        using namespace PHARE;
        using namespace std::literals;
        if (auto const e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::threaded::ScopeTimerMan::INSTANCE()
                .file_name(".phare_times.0.txt")
                // .force_strings()
                // .headers("fn"s, "dim"s, "layout"s, "alloc"s, "storage"s, "time"s)
                .init(); //
    })
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
    using enum LayoutMode;
    using Particles = typename Ions::particle_array_type;
    if constexpr (Particles::is_mapped)
        return construct_<IonUpdater<Ions, EM, GridLayout_t>>();
    else if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC))
        return construct_<IonUpdaterMultiPC<Ions, EM, GridLayout_t>>();
    else if constexpr (any_in(Particles::layout_mode, AoSTS, SoATS, SoAVXTS))
        return construct_<mkn::IonUpdaterMultiTS<Ions, EM, GridLayout_t>>();
    else
        return construct_<IonUpdaterPP<Ions, EM, GridLayout_t>>();
}


namespace detail::strings
{
    constexpr static std::string_view test    = "test,";
    constexpr static std::string_view update  = "update,";
    constexpr static std::string_view compare = "compare,";
    constexpr static std::string_view append  = "append,";
    constexpr static std::string_view convert = "convert,";
    constexpr static std::string_view cma     = ",";
} // namespace detail::strings


template<typename Patches>
void ref_update(UpdaterMode mode, Patches& patches)
{
    using Particles = Patches::value_type::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);
    for (auto& [layout, ions, em] : patches)
        get_updater_for(*ions, *em, layout).updatePopulations(*ions, *em, layout, dt, mode);
}

template<typename Patches>
void cmp_update(UpdaterMode mode, Patches& patches)
{
    using Particles = Patches::value_type::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);
    get_updater_for(*patches[0].ions, *patches[0].electromag, patches[0].layout)
        .updatePopulations(patches, dt, mode);
}


template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons_t<Particles_t, interp> ions{layout, "protons"};

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto seeding = [&]() {
        if constexpr (COMPARE_TO_ONE_PATCH)
            return seed;
        else
            return seeder++;
    }();

    auto disperse = [&](auto& particles) {
        delta_disperse(particles.domain_particles, seeding);
        delta_disperse(particles.patch_ghost_particles, seeding);
        if (shufle > 0)
        {
            shuffle(particles.domain_particles, seeding);
            shuffle(particles.patch_ghost_particles, seeding);
        }
    };

    auto add_particles = [&](auto& particles) {
        add_particles_in(particles.domain_particles, layout.AMRBox(), ppc);
        if constexpr (INCLUDE_PATCH_GHOST)
            add_ghost_particles(particles.patch_ghost_particles, layout.AMRBox(), ppc,
                                GridLayout_t::nbrParticleGhosts());
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);

    return ions;
}



template<typename Particles_t, typename GridLayout_t, typename Ions>
auto from_ions(GridLayout_t const& layout, Ions const& from)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons_t<Particles_t, interp> ions{layout, "protons"};
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto _add_particles_from = [&]<auto type>(auto& src, auto& dst) {
        // ParticleArrayService::reserve_ppc_in<type>(dst, ppc);
        append_particles<type>(src, dst /*, layout*/);
    };

    _add_particles_from.template operator()<ParticleType::Domain>(
        from.populations[0].particles.domain_particles,
        ions.populations[0].particles.domain_particles);
    if constexpr (INCLUDE_PATCH_GHOST)
        _add_particles_from.template operator()<ParticleType::Ghost>(
            from.populations[0].particles.patch_ghost_particles,
            ions.populations[0].particles.patch_ghost_particles);
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);
    return ions;
}


template<typename GridLayout_t, typename P0, typename P1>
void check_particles(GridLayout_t const& layout, P0& ref, P1& cmp_, std::size_t ghosts = 0)
{
    using CPU_ref
        = ParticleArray<P0::dimension,
                        ParticleArrayInternals<P0::dimension, LayoutMode::AoS, StorageMode::VECTOR,
                                               AllocatorMode::CPU, P0::impl>>;

    auto const box = grow(layout.AMRBox(), ghosts);
    auto const cmp = convert_particles_and_sort<CPU_ref>(cmp_, layout);
    sort_particles(ref, box);

    EXPECT_EQ(ref.size(), cmp.size());

    auto report = ref == cmp;
    if (report)
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays OK: " << P0::id() << " vs " << P1::id());
    }
    else
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays FAIL: " << P0::id() << " vs " << P1::id());
    }
    PHARE_LOG_LINE_STR("results: " << report.why());
    PHARE_LOG_LINE_STR("eg: " << ref[0]);

    EXPECT_EQ(ref, cmp);
}

template<bool GHOSTS = false, typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using ParticleArray_t = typename C::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::compare, ParticleArray_t::type_id,
                              detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    check_particles(layout, ref.populations[0].particles.domain_particles,
                    cmp.populations[0].particles.domain_particles);

    if constexpr (GHOSTS)
        check_particles(layout, ref.populations[0].particles.patch_ghost_particles,
                        cmp.populations[0].particles.patch_ghost_particles);

    using enum LayoutMode;
    using enum AllocatorMode;

    double diff = 1e-15;
    if constexpr (ParticleArray_t::alloc_mode == GPU_UNIFIED)
        diff *= 1e3; // atomics no order guaranteed#
    else if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS, SoATS))
        diff *= 1e1; // p2m op order diff

    {
        auto const freport
            = compare_tensor_fields(*ref.populations[0].F, *cmp.populations[0].F, diff);
        PHARE_LOG_LINE_STR("results: " << freport.why());
        EXPECT_TRUE(freport);
    }
    auto const rhoport = compare_fields(*ref.populations[0].rho, *cmp.populations[0].rho, diff);
    PHARE_LOG_LINE_STR("results: " << rhoport.why());
    EXPECT_TRUE(rhoport);
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

    using UsableElectromag_t = UsableElectromag<dim, alloc_mode>;
    using RefIons_t          = UsableIons_t<RefParticleArray_t, interp>;
    using CmpIons_t          = UsableIons_t<CmpParticleArray_t, interp>;

    GridLayout_t const layout{cells};


    MultiPatchIonUpdaterTest() {}

    void setup()
    {
        if constexpr (COMPARE_TO_ONE_PATCH)
            init_one();
        else
            init_n();

        ref_update(updater_mode, ref_patches);
    }

    void init_one()
    {
        auto& ref = ref_patches.emplace_back(layout, make_ions<RefParticleArray_t>(layout));
        for (std::size_t i = 0; i < n_patches; i++)
            cmp_patches.emplace_back(layout, from_ions<CmpParticleArray_t>(layout, ref.ions));
    }
    void init_n()
    {
        for (std::size_t i = 0; i < n_patches; i++)
        {
            auto& ref = ref_patches.emplace_back(layout, make_ions<RefParticleArray_t>(layout));
            cmp_patches.emplace_back(layout, from_ions<CmpParticleArray_t>(layout, ref.ions));
        }
    }

    void init_check() const
    {
        bool static constexpr init_value_check = 0;
        if constexpr (init_value_check)
        {
            auto ref = make_ions<RefParticleArray_t>(layout);
            auto cmp = from_ions<CmpParticleArray_t>(layout, *ref);
            compare(*this->layout, *ref, *cmp);
        }
    }

    void do_compare()
    {
        if constexpr (COMPARE_TO_ONE_PATCH)
            for (auto& cmp : cmp_patches)
                compare(*layout, ref_patches[0].ions, cmp.ions);
        else
            for (std::size_t i = 0; i < n_patches; i++)
                compare(*layout, ref_patches[i].ions, cmp_patches[i].ions);
    }


    template<typename Ions_t, typename ParticleArray_t_ = typename Ions_t::particle_array_type>
    struct Patch
    {
        using ParticleArray_t = ParticleArray_t_;
        using GridLayout_t    = MultiPatchIonUpdaterTest<Param>::GridLayout_t;
        using Electromag_t    = UsableElectromag_t::Super;

        GridLayout_t layout;
        Ions_t ions;
        UsableElectromag_t electromag{layout};
    };

    using RefPatch = Patch<RefIons_t>;
    using CmpPatch = Patch<CmpIons_t>;

    std::vector<aggregate_adapter<RefPatch>> ref_patches{};
    std::vector<aggregate_adapter<CmpPatch>> cmp_patches{};
};



// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     // TestParam<3, LayoutMode::SoAVXTS, AllocatorMode::CPU, 2, UpdaterMode::all>
     TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU, 2, UpdaterMode::all>
    // ,TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU, 2, UpdaterMode::all>
    // ,TestParam<3, LayoutMode::SoA, AllocatorMode::CPU, 2, UpdaterMode::all>

PHARE_WITH_MKN_GPU(

    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>
    // ,TestParam<3, LayoutMode::SoATS, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>

    // ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::domain_only>
    // ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::domain_only>
    // ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>
    // ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED, 2, UpdaterMode::all>

)

>;
// clang-format on

TYPED_TEST_SUITE(MultiPatchIonUpdaterTest, Permutations_t, );

template<typename MultiPatchIonUpdaterTest_t>
auto run(MultiPatchIonUpdaterTest_t& self)
{
    // using Particles = typename MultiPatchIonUpdaterTest_t::CmpParticleArray_t;
    // auto constexpr function_id
    //     = join_string_views_v<detail::strings::test, Particles::type_id, detail::strings::cma>;
    // PHARE_LOG_SCOPE(1, function_id);

    self.setup();

    cmp_update(MultiPatchIonUpdaterTest_t::updater_mode, self.cmp_patches);

    if (do_cmp)
        self.do_compare();
}




TYPED_TEST(MultiPatchIonUpdaterTest, updater_domain_only)
{
    run(*this);
}

template<std::size_t dim, typename internals> // used by gtest
void PrintTo(ParticleArray<dim, internals> const& arr, std::ostream* os)
{
    *os << arr;
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
