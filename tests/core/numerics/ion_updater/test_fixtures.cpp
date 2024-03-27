

#include "core/numerics/ion_updater/ion_updater.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"


struct IonUpdaterFixturesTest : public ::testing::Test
{
};

template<typename IonUpdater_t, typename Ions, typename EM, typename GridLayout_t>
auto construct_(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    PHARE::initializer::PHAREDict dict;
    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    return IonUpdater_t{dict["simulation"]["algo"]["ion_updater"]};
}

template<typename Ions, typename EM, typename GridLayout_t>
auto get_updater_for(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    using namespace PHARE::core;
    // there is only one, currently...
    return construct_<IonUpdater<Ions, EM, GridLayout_t>>(ions, em, layout);
}

template<typename Ions, typename EM, typename GridLayout_t>
void update(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    double dt = .001;
    get_updater_for(ions, em, layout)
        .updatePopulations(ions, em, layout, dt, PHARE::core::UpdaterMode::domain_only);
}

template<typename Particles_t, typename GridLayout_t>
auto evolve(GridLayout_t const& layout)
{
    using namespace PHARE::core;
    std::size_t constexpr static ppc = 2;

    auto constexpr static dim    = GridLayout_t::dimension;
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons_t<Particles_t, interp> ions{layout, "protons"};
    UsableElectromag<dim> em{layout};

    assert(ions.populations[0].particles.domain_particles.size() == 0);
    add_particles_in(ions.populations[0].particles.domain_particles, layout.AMRBox(), ppc);
    assert(ions.populations[0].particles.domain_particles.size() == layout.AMRBox().size() * ppc);

    update(*ions, *em, layout);

    assert(ions.populations[0].particles.domain_particles.size() > 0);

    return ions;
}


template<std::size_t dim, std::size_t cells = 5>
void test()
{
    using namespace PHARE::core;
    using PHARE_Types_t = PHARE_Types<dim, /*interp=*/1>;
    TestGridLayout<typename PHARE_Types_t::GridLayout_t> layout{cells};

    auto ref = evolve<PHARE::core::ParticleArray<dim>>(*layout);
}

TEST(IonUpdaterFixturesTest, works_1d)
{
    test<1>();
}

TEST(IonUpdaterFixturesTest, works_2d)
{
    test<2>();
}

TEST(IonUpdaterFixturesTest, works_3d)
{
    test<3>();
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}