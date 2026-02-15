//  tests/core/numerics/ion_updater/test_multi_updater.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

// USE HIP_VISIBLE_DEVICES OR CUDA_VISIBLE_DEVICES env vars

#define PHARE_UNDEF_ASSERT
#include "core/logger.hpp"
#define PHARE_SKIP_MPI_IN_CORE

#include "initializer/data_provider.hpp"
#include "core/utilities/thread_pool.hpp" // defaults to 1 thread, setup during static init!
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updaters.hpp" // IWYU pragma: keep
#include "core/data/particles/particle_array_appender.hpp"
#include "core/numerics/ion_updater/ion_updater_per_particle.hpp"

#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"


#include "gtest/gtest.h"
#include <cstddef>


namespace PHARE::core
{
// RUNTIME ENV VAR OVERRIDES
auto static const bytes     = get_env_as("PHARE_GPU_BYTES", std::uint64_t{500000000});
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{12});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{100});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1067});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const shufle    = get_env_as("PHARE_UNSORTED", std::size_t{0});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
auto static const n_threads = get_env_as("PHARE_THREADS", std::size_t{1});
auto static const cmp_only  = get_env_as("PHARE_CMP_ONLY", std::size_t{0});
auto static const ref_only  = get_env_as("PHARE_REF_ONLY", std::size_t{0});

bool static const premain = []() {
    assert(!(cmp_only and ref_only) && "Cant have both only");
    PHARE_WITH_MKN_GPU({ //
        ::mkn::gpu::setLimitMallocHeapSize(bytes);
    })

    PHARE_WITH_PHLOP({
        PHARE_LOG_LINE_STR("compare    : " << do_cmp);
        PHARE_LOG_LINE_STR("cells      : " << cells);
        PHARE_LOG_LINE_STR("ppc        : " << ppc);
        PHARE_LOG_LINE_STR("particles  : " << std::pow(cells, 3) * ppc);
        PHARE_LOG_LINE_STR("n_patches  : " << n_patches);
        PHARE_LOG_LINE_STR("seed       : " << seed);
        PHARE_LOG_LINE_SS("particle MB ≈ " << n_patches * std::pow(cells, 3) * ppc * 76 / 1e6);

        using namespace std::literals;
        if (auto const e = get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::threaded::ScopeTimerMan::INSTANCE().file_name(".phare_times.0.txt").init();
    })
    ThreadPool::threads_per_pool = n_threads;
    return true;
}();


auto& pool = *ThreadPool::INSTANCE().thread_pools[0];

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
    // else if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC))
    //     return construct_<IonUpdaterMultiPC<Ions, EM, GridLayout_t>>();
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
    using GridLayout_t = Patches::value_type::GridLayout_t;
    using Particles    = Patches::value_type::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);

    for (auto& [layout, ions, _, electromag] : patches)
    {
        auto updater       = get_updater_for(*ions, electromag, layout);
        using IonUpdater_t = std::decay_t<decltype(updater)>;
        using Boxing_t
            = std::decay_t<decltype(*selection_boxing_impl<IonUpdater_t, GridLayout_t>())>;

        Boxing_t const boxing{layout, {layout.AMRBox()}};
        updater.updatePopulations(*ions, electromag, boxing, dt, mode);
    }
}

template<typename Patches>
void cmp_update(UpdaterMode mode, Patches& patches)
{
    using GridLayout_t = Patches::value_type::GridLayout_t;
    using Boxing_t     = PHARE::core::UpdaterTileSetSelectionBoxing<GridLayout_t>;
    using Particles    = Patches::value_type::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);

    auto const& layout = patches[0].layout;
    std::unordered_map<std::string, Boxing_t> boxings;
    boxings.try_emplace("patch_id", Boxing_t{layout, {layout.AMRBox()}});
    get_updater_for(*patches[0].ions, patches[0].electromag, layout)
        .updatePopulations(patches, boxings, dt, mode);
}


template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons_t<Particles_t, interp> ions{layout, "protons"};

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto const disperse = [&](auto& particles) {
        delta_disperse(particles.domain_particles, seed);
        vary_velocity(particles.domain_particles, -6, 6, seed);
        if (shufle > 0)
            shuffle(particles.domain_particles, seed);


        delta_disperse(particles.level_ghost_particles, seed);
        vary_velocity(particles.level_ghost_particles, -6, 6, seed);
    };

    auto const particle_box = layout.AMRBox();

    auto add_particles = [&](auto& particles) {
        particles.domain_particles.reserve(particle_box.size() * ppc);
        add_particles_in(particles.domain_particles, particle_box, ppc);

        add_ghost_particles(particles.level_ghost_particles, particle_box, ppc,
                            GridLayout_t::nbrParticleGhosts());
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);


    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);

    return ions;
}



template<typename Particles_t, typename GridLayout_t, typename Ions>
auto from_ions(GridLayout_t const& layout, Ions const& from)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons_t<Particles_t, interp> ions{layout, "protons"};
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto _add_particles_from = [&]<auto type>(auto& src, auto& dst) {
        ParticleArrayService::reserve_ppc_in<type>(dst, ppc);
        append_particles<type>(src, dst /*, layout*/);
    };

    if constexpr (Particles_t::layout_mode == AoSTS)
    {
        PHARE_LOG_LINE_SS("n tiles: " << ions.populations[0].particles.domain_particles().size());
    }

    _add_particles_from.template operator()<ParticleType::Domain>(
        from.populations[0].particles.domain_particles,
        ions.populations[0].particles.domain_particles);

    _add_particles_from.template operator()<ParticleType::LevelGhost>(
        from.populations[0].particles.level_ghost_particles,
        ions.populations[0].particles.level_ghost_particles);

    auto const particle_box = layout.AMRBox();
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);
    return std::move(ions);
}


template<typename GridLayout_t, typename P0, typename P1>
void check_particles(GridLayout_t const& layout, P0& ref, P1& cmp_, double const atol)
{
    using CPU_ref = ParticleArray<ParticleArrayOptions{P0::dimension, LayoutMode::AoS,
                                                       StorageMode::VECTOR, AllocatorMode::CPU}>;

    auto const box = layout.AMRBox();
    particle_array_domain_is_valid(cmp_, box);

    auto const cmp = convert_particles_and_sort<CPU_ref>(cmp_, layout);
    sort_particles(ref, box);

    EXPECT_EQ(ref.size(), cmp.size());

    auto const report = compare_particles(ref, cmp, atol);
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

    // EXPECT_EQ(ref, cmp);
}

template<typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using ParticleArray_t = C::ParticleArray_t;
    auto constexpr function_id
        = join_string_views_v<detail::strings::compare, ParticleArray_t::type_id,
                              detail::strings::cma>;
    PHARE_LOG_SCOPE(1, function_id);

    using enum LayoutMode;
    using enum AllocatorMode;

    double diff = 1e-15;
    if constexpr (ParticleArray_t::alloc_mode == GPU_UNIFIED)
        diff *= 1e3; // atomics no order guaranteed
    else if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS, SoATS))
        diff *= 1e1; // p2m op order diff

    check_particles(layout, ref.populations[0].particles.domain_particles,
                    cmp.populations[0].particles.domain_particles, diff);

    zero_ghost_layer(layout, ref.populations[0].F, cmp.populations[0].F);
    auto const freport
        = compare_reduced_tensor_fields(ref.populations[0].F, cmp.populations[0].F, diff);
    PHARE_LOG_LINE_STR("results: " << freport.why());
    EXPECT_TRUE(freport);

    zero_ghost_layer(layout, ref.populations[0].particleDensity_,
                     cmp.populations[0].particleDensity_);
    auto const rhoPport = compare_reduced_fields(ref.populations[0].particleDensity_,
                                                 cmp.populations[0].particleDensity_, diff);
    PHARE_LOG_LINE_STR("results: " << rhoPport.why());
    EXPECT_TRUE(rhoPport);

    zero_ghost_layer(layout, ref.populations[0].chargeDensity_, cmp.populations[0].chargeDensity_);
    auto const rhoCport = compare_reduced_fields(ref.populations[0].chargeDensity_,
                                                 cmp.populations[0].chargeDensity_, diff);
    PHARE_LOG_LINE_STR("results: " << rhoCport.why());
    EXPECT_TRUE(rhoCport);
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode, auto _updater_mode>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), PHARE::AllocatorMode>);
    auto constexpr static dim          = _dim;
    auto constexpr static layout_mode  = _layout_mode;
    auto constexpr static alloc_mode   = _alloc_mode;
    auto constexpr static updater_mode = _updater_mode;
};



template<typename Param>
struct MultiPatchIonUpdaterTest : public ::testing::Test
{
    auto constexpr static dim          = Param::dim;
    auto constexpr static interp       = 1;
    auto constexpr static layout_mode  = Param::layout_mode;
    auto constexpr static alloc_mode   = Param::alloc_mode;
    auto constexpr static updater_mode = Param::updater_mode;


    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::GridLayout_t>;
    using RefParticleArray_t = AoSMappedParticleArray<dim>;
    using CmpParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    using UsableElectromag_t
        = UsableElectromag<typename GridLayout_t::Super, alloc_mode, layout_mode>;
    using RefIons_t = UsableIons_t<RefParticleArray_t, interp>;
    using RefEM_t   = UsableElectromag<typename GridLayout_t::Super>;
    using CmpIons_t = UsableIons_t<CmpParticleArray_t, interp>;
    using CmpEM_t   = UsableElectromag<typename GridLayout_t::Super, alloc_mode, layout_mode>;

    GridLayout_t const layout{cells};


    MultiPatchIonUpdaterTest() {}

    void run()
    {
        cmp_patches.reserve(n_patches);
        auto& ref = ref_patches.emplace_back(layout, make_ions<RefParticleArray_t>(layout));

        if (!ref_only)
            for (std::size_t i = 0; i < n_patches; i++)
                cmp_patches.emplace_back(layout, from_ions<CmpParticleArray_t>(layout, ref.ions));

        if (!cmp_only)
            ref_update(updater_mode, ref_patches);

        if (!ref_only)
            cmp_update(updater_mode, cmp_patches);

        if (do_cmp)
            for (auto& cmp : cmp_patches)
                compare(*layout, ref_patches[0].ions, cmp.ions);
    }


    template<typename Ions_t, typename EM>
    struct Patch
    {
        using ParticleArray_t = Ions_t::particle_array_type;
        using GridLayout_t    = MultiPatchIonUpdaterTest<Param>::GridLayout_t;
        using Electromag_t    = EM::Super;

        GridLayout_t layout;
        Ions_t ions;
        EM em{layout};
        Electromag_t electromag = *em;

        std::string patchID() const { return "patch_id"; }
    };


    using RefPatch = Patch<RefIons_t, RefEM_t>;
    using CmpPatch = Patch<CmpIons_t, CmpEM_t>;

    std::vector<aggregate_adapter<RefPatch>> ref_patches{};
    std::vector<aggregate_adapter<CmpPatch>> cmp_patches{};
};



// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU, UpdaterMode::domain_only>
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU, UpdaterMode::all>


PHARE_WITH_GPU(

    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, UpdaterMode::domain_only>
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, UpdaterMode::all>

)

>;
// clang-format on

TYPED_TEST_SUITE(MultiPatchIonUpdaterTest, Permutations_t, );

TYPED_TEST(MultiPatchIonUpdaterTest, updater_domain_only)
{
    this->run();
}

template<auto opts, typename internals> // used by gtest
void PrintTo(ParticleArray<opts, internals> const& arr, std::ostream* os)
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
    PHARE::core::MemoryMonitoring::PRINT();
    return r;
}
