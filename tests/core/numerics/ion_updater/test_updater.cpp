
#include "phare_core.hpp"

#include "core/utilities/box/box.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"


#include "gtest/gtest.h"

using namespace PHARE::core;



using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return density(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

Return vx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vy(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vthx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vthy(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vthz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return bx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return by(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return bz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}




std::size_t nbrPartPerCell = 1000;

using InitFunctionT = PHARE::initializer::InitFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};

    dict["ions"]["nbrPopulations"]                          = std::size_t{2};
    dict["ions"]["pop0"]["name"]                            = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                            = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"]
        = static_cast<int>(nbrPartPerCell);
    dict["ions"]["pop0"]["particle_initializer"]["charge"] = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]  = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]                            = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]                            = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"]
        = static_cast<int>(nbrPartPerCell);
    dict["ions"]["pop1"]["particle_initializer"]["charge"] = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]  = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    return dict;
}
static auto init_dict = createDict();



template<std::size_t dim, std::size_t interporder>
struct DimInterp
{
    static constexpr auto dimension    = dim;
    static constexpr auto interp_order = interporder;
};



// the Electromag and Ions used in this test
// need their resources pointers (Fields and ParticleArrays) to set manually
// to buffers. ElectromagBuffer and IonsBuffer encapsulate these buffers



template<std::size_t dim, std::size_t interp_order>
struct ElectromagBuffers
{
    using PHARETypes       = PHARE::core::PHARE_Types<dim, interp_order>;
    using Grid             = typename PHARETypes::Grid_t;
    using GridLayout       = typename PHARETypes::GridLayout_t;
    using Electromag       = typename PHARETypes::Electromag_t;
    using UsableVecFieldND = UsableVecField<dim>;

    UsableVecFieldND B, E;

    ElectromagBuffers(GridLayout const& layout)
        : B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
    {
    }

    ElectromagBuffers(ElectromagBuffers const& source, GridLayout const& layout)
        : ElectromagBuffers{layout}
    {
        B.copyData(source.B);
        E.copyData(source.E);
    }


    void setBuffers(Electromag& EM)
    {
        B.set_on(EM.B);
        E.set_on(EM.E);
    }
};




template<std::size_t dim, std::size_t interp_order>
struct IonsBuffers
{
    using PHARETypes                 = PHARE::core::PHARE_Types<dim, interp_order>;
    using UsableVecFieldND           = UsableVecField<dim>;
    using Grid                       = typename PHARETypes::Grid_t;
    using GridLayout                 = typename PHARETypes::GridLayout_t;
    using Ions                       = typename PHARETypes::Ions_t;
    using ParticleArray              = typename PHARETypes::ParticleArray_t;
    using ParticleInitializerFactory = typename PHARETypes::ParticleInitializerFactory;

    Grid ionDensity;
    Grid ionMassDensity;
    Grid protonDensity;
    Grid alphaDensity;

    UsableVecFieldND protonF, alphaF, Vi;
    UsableTensorField<dim> M, alpha_M, protons_M;

    static constexpr int ghostSafeMapLayer = ghostWidthForParticles<interp_order>() + 1;

    ParticleArray protonDomain;
    ParticleArray protonPatchGhost;
    ParticleArray protonLevelGhost;
    ParticleArray protonLevelGhostOld;
    ParticleArray protonLevelGhostNew;

    ParticleArray alphaDomain;
    ParticleArray alphaPatchGhost;
    ParticleArray alphaLevelGhost;
    ParticleArray alphaLevelGhostOld;
    ParticleArray alphaLevelGhostNew;

    ParticlesPack<ParticleArray> protonPack;
    ParticlesPack<ParticleArray> alphaPack;

    IonsBuffers(GridLayout const& layout)
        : ionDensity{"rho", HybridQuantity::Scalar::rho,
                     layout.allocSize(HybridQuantity::Scalar::rho)}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonDensity{"protons_rho", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaDensity{"alpha_rho", HybridQuantity::Scalar::rho,
                       layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonF{"protons_flux", layout, HybridQuantity::Vector::V}
        , alphaF{"alpha_flux", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , alpha_M{"alpha_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protons_M{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protonDomain{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonPatchGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhostOld{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhostNew{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaDomain{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaPatchGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhostOld{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhostNew{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonPack{"protons",         &protonDomain,        &protonPatchGhost,
                     &protonLevelGhost, &protonLevelGhostOld, &protonLevelGhostNew}
        , alphaPack{"alpha",          &alphaDomain,        &alphaPatchGhost,
                    &alphaLevelGhost, &alphaLevelGhostOld, &alphaLevelGhostNew}
    {
    }


    IonsBuffers(IonsBuffers const& source, GridLayout const& layout)
        : ionDensity{"rho", HybridQuantity::Scalar::rho,
                     layout.allocSize(HybridQuantity::Scalar::rho)}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonDensity{"protons_rho", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaDensity{"alpha_rho", HybridQuantity::Scalar::rho,
                       layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonF{"protons_flux", layout, HybridQuantity::Vector::V}
        , alphaF{"alpha_flux", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , alpha_M{"alpha_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protons_M{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protonDomain{source.protonDomain}
        , protonPatchGhost{source.protonPatchGhost}
        , protonLevelGhost{source.protonLevelGhost}
        , protonLevelGhostOld{source.protonLevelGhostOld}
        , protonLevelGhostNew{source.protonLevelGhostNew}
        , alphaDomain{source.alphaDomain}
        , alphaPatchGhost{source.alphaPatchGhost}
        , alphaLevelGhost{source.alphaLevelGhost}
        , alphaLevelGhostOld{source.alphaLevelGhostOld}
        , alphaLevelGhostNew{source.alphaLevelGhostNew}
        , protonPack{"protons",         &protonDomain,        &protonPatchGhost,
                     &protonLevelGhost, &protonLevelGhostOld, &protonLevelGhostNew}
        , alphaPack{"alpha",          &alphaDomain,        &alphaPatchGhost,
                    &alphaLevelGhost, &alphaLevelGhostOld, &alphaLevelGhostNew}

    {
        ionDensity.copyData(source.ionDensity);
        ionMassDensity.copyData(source.ionMassDensity);
        protonDensity.copyData(source.protonDensity);
        alphaDensity.copyData(source.alphaDensity);

        protonF.copyData(source.protonF);
        alphaF.copyData(source.alphaF);
        Vi.copyData(source.Vi);
    }

    void setBuffers(Ions& ions)
    {
        {
            auto const& [V, m, d, md] = ions.getCompileTimeResourcesViewList();
            Vi.set_on(V);
            M.set_on(m);
            d.setBuffer(&ionDensity);
            md.setBuffer(&ionMassDensity);
        }

        auto& pops = ions.getRunTimeResourcesViewList();
        {
            auto const& [F, M, d, particles] = pops[0].getCompileTimeResourcesViewList();
            d.setBuffer(&protonDensity);
            protons_M.set_on(M);
            protonF.set_on(F);
            particles.setBuffer(&protonPack);
        }

        {
            auto const& [F, M, d, particles] = pops[1].getCompileTimeResourcesViewList();
            d.setBuffer(&alphaDensity);
            alpha_M.set_on(M);
            alphaF.set_on(F);
            particles.setBuffer(&alphaPack);
        }
    }
};




template<typename DimInterpT>
struct IonUpdaterTest : public ::testing::Test
{
    static constexpr auto dim          = DimInterpT::dimension;
    static constexpr auto interp_order = DimInterpT::interp_order;
    using PHARETypes                   = PHARE::core::PHARE_Types<dim, interp_order>;
    using Ions                         = typename PHARETypes::Ions_t;
    using Electromag                   = typename PHARETypes::Electromag_t;
    using GridLayout    = typename PHARE::core::GridLayout<GridLayoutImplYee<dim, interp_order>>;
    using ParticleArray = typename PHARETypes::ParticleArray_t;
    using ParticleInitializerFactory = typename PHARETypes::ParticleInitializerFactory;

    using IonUpdater = typename PHARE::core::IonUpdater<Ions, Electromag, GridLayout>;
    using Boxing_t   = PHARE::core::UpdaterSelectionBoxing<IonUpdater, GridLayout>;


    double dt{0.01};

    // grid configuration
    std::array<int, dim> ncells;
    GridLayout layout;
    // assumes no level ghost cells

    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout::nbrParticleGhosts())}};


    // data for electromagnetic fields
    using Field            = typename PHARETypes::Grid_t;
    using VecField         = typename PHARETypes::VecField_t;
    using UsableVecFieldND = UsableVecField<dim>;

    ElectromagBuffers<dim, interp_order> emBuffers;
    IonsBuffers<dim, interp_order> ionsBuffers;

    Electromag EM{init_dict["electromag"]};
    Ions ions{init_dict["ions"]};



    IonUpdaterTest()
        : ncells{100}
        , layout{{0.1}, {100u}, {{0.}}}
        , emBuffers{layout}
        , ionsBuffers{layout}
    {
        emBuffers.setBuffers(EM);
        ionsBuffers.setBuffers(ions);


        // ok all resources pointers are set to buffers
        // now let's initialize Electromag fields to user input functions
        // and ion population particles to user supplied moments


        EM.initialize(layout);
        for (auto& pop : ions)
        {
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
        }


        // now all domain particles are loaded we need to manually insert
        // ghost particles (this is in reality SAMRAI's job)
        // these are needed if we want all used nodes to be complete


        // in 1D we assume left border is touching the level border
        // and right is touching another patch
        // so on the left no patchGhost but levelGhost(and old and new)
        // on the right no levelGhost but patchGhosts


        for (auto& pop : ions)
        {
            if constexpr (dim == 1)
            {
                int firstPhysCell = layout.physicalStartIndex(QtyCentering::dual, Direction::X);
                int lastPhysCell  = layout.physicalEndIndex(QtyCentering::dual, Direction::X);
                auto firstAMRCell = layout.localToAMR(Point{firstPhysCell});
                auto lastAMRCell  = layout.localToAMR(Point{lastPhysCell});

                // we need to put levelGhost particles in the cell just to the
                // left of the first cell. In reality these particles should
                // come from splitting particles of the next coarser level
                // in this test we just copy those of the first cell
                // we also assume levelGhostOld and New are the same particles
                // for simplicity

                auto& domainPart        = pop.domainParticles();
                auto& levelGhostPartOld = pop.levelGhostParticlesOld();
                auto& levelGhostPartNew = pop.levelGhostParticlesNew();
                auto& levelGhostPart    = pop.levelGhostParticles();
                auto& patchGhostPart    = pop.patchGhostParticles();


                // copies need to be put in the ghost cell
                // we have copied particles be now their iCell needs to be udpated
                // our choice is :
                //
                // first order:
                //
                //   ghost| domain...
                // [-----]|[-----][-----][-----][-----][-----]
                //     ^      v
                //     |      |
                //     -------|
                //
                // second and third order:

                //   ghost        | domain...
                // [-----]|[-----][-----][-----][-----][-----][-----]
                //     ^      ^       v     v
                //     |      |       |     |
                //     -------|-------|     |
                //            ---------------
                for (auto const& part : domainPart)
                {
                    if constexpr (interp_order == 2 or interp_order == 3)
                    {
                        if (part.iCell[0] == firstAMRCell[0]
                            or part.iCell[0] == firstAMRCell[0] + 1)
                        {
                            auto p{part};
                            p.iCell[0] -= 2;
                            levelGhostPartOld.push_back(p);
                        }
                    }
                    else if constexpr (interp_order == 1)
                    {
                        if (part.iCell[0] == firstAMRCell[0])
                        {
                            auto p{part};
                            p.iCell[0] -= 1;
                            levelGhostPartOld.push_back(p);
                        }
                    }
                }


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPartNew));


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPart));

                std::size_t const ghosts_cells = (interp_order == 1 ? 1 : 2);
                EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
                EXPECT_EQ(levelGhostPartOld.size(), ghosts_cells * nbrPartPerCell);
                EXPECT_EQ(patchGhostPart.size(), 0);


            } // end 1D
        } // end pop loop
        PHARE::core::depositParticles(ions, layout, Interpolator<dim, interp_order>{},
                                      PHARE::core::DomainDeposit{});


        PHARE::core::depositParticles(ions, layout, Interpolator<dim, interp_order>{},
                                      PHARE::core::LevelGhostDeposit{});


        ions.computeDensity();
        ions.computeBulkVelocity();
    } // end Ctor



    void fillIonsMomentsGhosts()
    {
        using Interpolator = typename IonUpdater::Interpolator;
        Interpolator interpolate;

        for (auto& pop : this->ions)
        {
            double alpha = 0.5;
            interpolate(makeIndexRange(pop.levelGhostParticlesNew()), pop.density(), pop.flux(),
                        layout,
                        /*coef = */ alpha);


            interpolate(makeIndexRange(pop.levelGhostParticlesOld()), pop.density(), pop.flux(),
                        layout,
                        /*coef = */ (1. - alpha));
        }
    }



    void checkMomentsHaveEvolved(IonsBuffers<dim, interp_order> const& ionsBufferCpy)
    {
        auto& populations = this->ions.getRunTimeResourcesViewList();

        auto& protonDensity = populations[0].density();
        auto& protonFx      = populations[0].flux().getComponent(Component::X);
        auto& protonFy      = populations[0].flux().getComponent(Component::Y);
        auto& protonFz      = populations[0].flux().getComponent(Component::Z);

        auto& alphaDensity = populations[1].density();
        auto& alphaFx      = populations[1].flux().getComponent(Component::X);
        auto& alphaFy      = populations[1].flux().getComponent(Component::Y);
        auto& alphaFz      = populations[1].flux().getComponent(Component::Z);

        auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

        auto nonZero = [&](auto const& field) {
            auto sum = 0.;
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                sum += std::abs(field(ix));
            }
            EXPECT_GT(sum, 0.);
        };

        auto check = [&](auto const& newField, auto const& originalField) {
            nonZero(newField);
            nonZero(originalField);
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                auto evolution = std::abs(newField(ix) - originalField(ix));
                //  should check that moments are still compatible with user inputs also
                EXPECT_TRUE(evolution > 0.0);
                if (evolution <= 0.0)
                    std::cout << "after update : " << newField(ix)
                              << " before update : " << originalField(ix)
                              << " evolution : " << evolution << " ix : " << ix << "\n";
            }
        };

        check(protonDensity, ionsBufferCpy.protonDensity);


        check(protonFx, ionsBufferCpy.protonF(Component::X));
        check(protonFy, ionsBufferCpy.protonF(Component::Y));
        check(protonFz, ionsBufferCpy.protonF(Component::Z));

        check(alphaDensity, ionsBufferCpy.alphaDensity);
        check(alphaFx, ionsBufferCpy.alphaF(Component::X));
        check(alphaFy, ionsBufferCpy.alphaF(Component::Y));
        check(alphaFz, ionsBufferCpy.alphaF(Component::Z));

        check(ions.density(), ionsBufferCpy.ionDensity);
        check(ions.velocity().getComponent(Component::X), ionsBufferCpy.Vi(Component::X));
        check(ions.velocity().getComponent(Component::Y), ionsBufferCpy.Vi(Component::Y));
        check(ions.velocity().getComponent(Component::Z), ionsBufferCpy.Vi(Component::Z));
    }



    void checkDensityIsAsPrescribed()
    {
        auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

        auto check = [&](auto const& density, auto const& function) {
            std::vector<std::size_t> ixes;
            std::vector<double> x;

            for (auto ix = ix0; ix < ix1; ++ix)
            {
                ixes.emplace_back(ix);
                x.emplace_back(layout.cellCenteredCoordinates(ix)[0]);
            }

            auto functionXPtr = function(x); // keep alive
            EXPECT_EQ(functionXPtr->size(), (ix1 - ix0));

            auto& functionX = *functionXPtr;

            for (std::size_t i = 0; i < functionX.size(); i++)
            {
                auto ix   = ixes[i];
                auto diff = std::abs(density(ix) - functionX[i]);

                EXPECT_GE(0.07, diff);

                if (diff >= 0.07)
                    std::cout << "actual : " << density(ix) << " prescribed : " << functionX[i]
                              << " diff : " << diff << " ix : " << ix << "\n";
            }
        };

        auto& populations   = this->ions.getRunTimeResourcesViewList();
        auto& protonDensity = populations[0].density();
        auto& alphaDensity  = populations[1].density();

        check(protonDensity, density);
        check(alphaDensity, density);
    }
};



using DimInterps = ::testing::Types<DimInterp<1, 1>, DimInterp<1, 2>, DimInterp<1, 3>>;


TYPED_TEST_SUITE(IonUpdaterTest, DimInterps, );


#if 0

TYPED_TEST(IonUpdaterTest, ionUpdaterTakesPusherParamsFromPHAREDictAtConstruction)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};
}


// the following 3 tests are testing the fixture is well configured.


TYPED_TEST(IonUpdaterTest, loadsDomainPatchAndLevelGhostParticles)
{
    auto check = [this](std::size_t nbrGhostCells, auto& pop) {
        EXPECT_EQ(this->layout.nbrCells()[0] * nbrPartPerCell, pop.domainParticles().size());
        EXPECT_EQ(0, pop.patchGhostParticles().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticlesOld().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticlesNew().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticles().size());
    };


    if constexpr (TypeParam::dimension == 1)
    {
        for (auto& pop : this->ions)
        {
            if constexpr (TypeParam::interp_order == 1)
            {
                check(1, pop);
            }
            else if constexpr (TypeParam::interp_order == 2 or TypeParam::interp_order == 3)
            {
                check(2, pop);
            }
        }
    }
}




TYPED_TEST(IonUpdaterTest, loadsLevelGhostParticlesOnLeftGhostArea)
{
    int firstPhysCell = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto firstAMRCell = this->layout.localToAMR(Point{firstPhysCell});

    if constexpr (TypeParam::dimension == 1)
    {
        for (auto& pop : this->ions)
        {
            if constexpr (TypeParam::interp_order == 1)
            {
                for (auto const& part : pop.levelGhostParticles())
                {
                    EXPECT_EQ(firstAMRCell[0] - 1, part.iCell[0]);
                }
            }
            else if constexpr (TypeParam::interp_order == 2 or TypeParam::interp_order == 3)
            {
                typename IonUpdaterTest<TypeParam>::ParticleArray copy{pop.levelGhostParticles()};
                auto firstInOuterMostCell = std::partition(
                    std::begin(copy), std::end(copy), [&firstAMRCell](auto const& particle) {
                        return particle.iCell[0] == firstAMRCell[0] - 1;
                    });
                EXPECT_EQ(nbrPartPerCell, std::distance(std::begin(copy), firstInOuterMostCell));
                EXPECT_EQ(nbrPartPerCell, std::distance(firstInOuterMostCell, std::end(copy)));
            }
        }
    }
}




// start of PHARE TESTS



TYPED_TEST(IonUpdaterTest, particlesUntouchedInMomentOnlyMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
                                 UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);


    auto& populations = this->ions.getRunTimeResourcesViewList();

    auto checkIsUnTouched = [](auto const& original, auto const& cpy) {
        // no particles should have moved, so none should have left the domain
        EXPECT_EQ(cpy.size(), original.size());
        for (std::size_t iPart = 0; iPart < original.size(); ++iPart)
        {
            EXPECT_EQ(cpy[iPart].iCell[0], original[iPart].iCell[0]);
            EXPECT_DOUBLE_EQ(cpy[iPart].delta[0], original[iPart].delta[0]);

            for (std::size_t iDir = 0; iDir < 3; ++iDir)
            {
                EXPECT_DOUBLE_EQ(cpy[iPart].v[iDir], original[iPart].v[iDir]);
            }
        }
    };

    checkIsUnTouched(populations[0].patchGhostParticles(), ionsBufferCpy.protonPatchGhost);
    checkIsUnTouched(populations[0].levelGhostParticles(), ionsBufferCpy.protonLevelGhost);
    checkIsUnTouched(populations[0].levelGhostParticlesOld(), ionsBufferCpy.protonLevelGhostOld);
    checkIsUnTouched(populations[0].levelGhostParticlesNew(), ionsBufferCpy.protonLevelGhostNew);

    checkIsUnTouched(populations[1].patchGhostParticles(), ionsBufferCpy.alphaPatchGhost);
    checkIsUnTouched(populations[1].levelGhostParticles(), ionsBufferCpy.alphaLevelGhost);
    checkIsUnTouched(populations[1].levelGhostParticlesOld(), ionsBufferCpy.alphaLevelGhost);
    checkIsUnTouched(populations[1].levelGhostParticlesNew(), ionsBufferCpy.alphaLevelGhost);
}




// TYPED_TEST(IonUpdaterTest, particlesAreChangedInParticlesAndMomentsMode)
//{
//    typename IonUpdaterTest<TypeParam>::IonUpdater
//    ionUpdater{init_dict["simulation"]["pusher"]};
//
//    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};
//
//    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
//                                 UpdaterMode::particles_and_moments);
//
//    this->fillIonsMomentsGhosts();
//
//    ionUpdater.updateIons(this->ions);
//
//    auto& populations = this->ions.getRunTimeResourcesViewList();
//
//    EXPECT_NE(ionsBufferCpy.protonDomain.size(), populations[0].domainParticles().size());
//    EXPECT_NE(ionsBufferCpy.alphaDomain.size(), populations[1].domainParticles().size());
//
//    // cannot think of anything else to check than checking that the number of particles
//    // in the domain have changed after them having been pushed.
//}



TYPED_TEST(IonUpdaterTest, momentsAreChangedInParticlesAndMomentsMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt, UpdaterMode::all);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}




TYPED_TEST(IonUpdaterTest, momentsAreChangedInMomentsOnlyMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
                                 UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}



TYPED_TEST(IonUpdaterTest, thatNoNaNsExistOnPhysicalNodesMoments)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
                                 UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto& pop : this->ions)
    {
        for (auto ix = ix0; ix <= ix1; ++ix)
        {
            auto& density = pop.density();
            auto& flux    = pop.flux();

            auto& fx = flux.getComponent(Component::X);
            auto& fy = flux.getComponent(Component::Y);
            auto& fz = flux.getComponent(Component::Z);

            EXPECT_FALSE(std::isnan(density(ix)));
            EXPECT_FALSE(std::isnan(fx(ix)));
            EXPECT_FALSE(std::isnan(fy(ix)));
            EXPECT_FALSE(std::isnan(fz(ix)));
        }
    }
}
#endif


template<std::size_t dim>
void pseudo_randomize_deltas(auto& parts)
{
    std::size_t k = 0;
    for (auto& p : parts)
        for (std::uint16_t i = 0; i < dim; ++i)
        {
            if (k == 9)
                k = 0;
            double dir = k % 2 == 0 ? -1 : 1;
            p.delta[i] += dir * .05 * k;
            assert(p.delta[i] > 0 and p.delta[i] < 1);
            ++k;
        }
}

template<typename GridLayout>
void populate_particles(auto& ions, GridLayout& layout)
{
    auto constexpr dim = GridLayout::dimension;
    Particle<dim> const particle{1.0 / nbrPartPerCell, 1, ConstArray<int, dim>(),
                                 ConstArray<double, dim>(.5), ConstArray<double, 3>(1)};

    auto const shift_particle = [](auto p, auto const bix) {
        p.iCell = bix;
        return p;
    };

    auto const ghost_box = grow(layout.AMRBox(), 1);
    for (auto& pop : ions)
    {
        pop.domainParticles().clear();
        pop.levelGhostParticles().clear();
        pop.patchGhostParticles().clear();

        EXPECT_TRUE(pop.domainParticles().is_consistent());
        EXPECT_TRUE(pop.levelGhostParticles().is_consistent());

        for (auto const& bix : layout.AMRBox())
            for (std::size_t i = 0; i < nbrPartPerCell; ++i)
                pop.domainParticles().emplace_back(shift_particle(particle, *bix));

        for (auto const& bix : ghost_box)
            if (not isIn(bix, layout.AMRBox()))
                for (std::size_t i = 0; i < nbrPartPerCell; ++i)
                    pop.levelGhostParticles().emplace_back(shift_particle(particle, *bix));

        for (auto const& p : pop.levelGhostParticles())
        {
            EXPECT_TRUE(not isIn(p, layout.AMRBox()));
        }

        pseudo_randomize_deltas<dim>(pop.domainParticles());
        pseudo_randomize_deltas<dim>(pop.levelGhostParticles());

        EXPECT_TRUE(pop.domainParticles().is_consistent());
        EXPECT_TRUE(pop.levelGhostParticles().is_consistent());
    }
}


TYPED_TEST(IonUpdaterTest, hardCodedRegressionTest0)
{
    using IonUpdater      = TestFixture::IonUpdater;
    auto constexpr dim    = TestFixture::dim;
    auto constexpr interp = TestFixture::interp_order;

    std::array<double, 3> expectations
        = {199.99999999999889, 199.99999999999659, 199.99999999999673};


    populate_particles(this->ions, this->layout);


    for (auto& pop : this->ions)
    {
        EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
        EXPECT_EQ(pop.levelGhostParticles().size(), 2 * nbrPartPerCell);
        EXPECT_EQ(pop.patchGhostParticles().size(), 0);
    }

    IonUpdater ionUpdater{init_dict["simulation"]["algo"]["ion_updater"]};
    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
                                 UpdaterMode::domain_only);

    std::int64_t icellSum = 0;
    double densitySum = 0, fluxXSum = 0, fluxYSum = 0, fluxZSum = 0;
    for (auto& pop : this->ions)
    {
        densitySum += PHARE::core::sum(pop.density());
        fluxXSum += PHARE::core::sum(pop.flux()[0]);
        fluxYSum += PHARE::core::sum(pop.flux()[1]);
        fluxZSum += PHARE::core::sum(pop.flux()[2]);
        icellSum += PHARE::core::sum_from(pop.domainParticles(),
                                          [](auto const& p) { return PHARE::core::sum(p.iCell); });
    }


    EXPECT_EQ(icellSum, 9900000);
    EXPECT_DOUBLE_EQ(densitySum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxXSum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxYSum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxZSum, expectations[interp - 1]);
}



TYPED_TEST(IonUpdaterTest, hardCodedRegressionTest1)
{
    using IonUpdater      = TestFixture::IonUpdater;
    auto constexpr dim    = TestFixture::dim;
    auto constexpr interp = TestFixture::interp_order;

    std::array<double, 3> expectations
        = {199.99999999999889, 199.99999999999659, 199.99999999999673};

    populate_particles(this->ions, this->layout);

    double deltaSum = 0;
    for (auto& pop : this->ions)
    {
        EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
        EXPECT_EQ(pop.levelGhostParticles().size(), 2 * nbrPartPerCell);
        EXPECT_EQ(pop.patchGhostParticles().size(), 0);
        deltaSum += PHARE::core::sum_from(pop.domainParticles(),
                                          [](auto const& p) { return PHARE::core::sum(p.delta); });
    }
    EXPECT_DOUBLE_EQ(deltaSum, 103333.29999993632);

    IonUpdater ionUpdater{init_dict["simulation"]["algo"]["ion_updater"]};
    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt, UpdaterMode::all);

    for (auto& pop : this->ions)
    {
        EXPECT_EQ(pop.domainParticles().size(), 100000);
        std::size_t domainInDomain = 0;
        for (auto const& p : pop.domainParticles())
            domainInDomain += isIn(p, this->layout.AMRBox()) ? 1 : 0;
        EXPECT_EQ(domainInDomain, pop.domainParticles().size());
        EXPECT_EQ(pop.levelGhostParticles().size(), 2 * nbrPartPerCell);
        EXPECT_EQ(pop.patchGhostParticles().size(), 0);
    }

    std::size_t icellSum = 0, pCount = 0;
    deltaSum          = 0;
    double densitySum = 0, fluxXSum = 0, fluxYSum = 0, fluxZSum = 0;
    for (auto& pop : this->ions)
    {
        densitySum += PHARE::core::sum(pop.density());
        fluxXSum += PHARE::core::sum(pop.flux()[0]);
        fluxYSum += PHARE::core::sum(pop.flux()[1]);
        fluxZSum += PHARE::core::sum(pop.flux()[2]);
        icellSum += PHARE::core::sum_from(pop.domainParticles(),
                                          [](auto const& p) { return PHARE::core::sum(p.iCell); });
        deltaSum += PHARE::core::sum_from(pop.domainParticles(),
                                          [](auto const& p) { return PHARE::core::sum(p.delta); });
        pCount += pop.domainParticles().size();
    }

    EXPECT_EQ(pCount, 200000);
    EXPECT_EQ(icellSum, 9900000);
    EXPECT_DOUBLE_EQ(deltaSum, 123333.30000002828);
    EXPECT_DOUBLE_EQ(densitySum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxXSum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxYSum, expectations[interp - 1]);
    EXPECT_DOUBLE_EQ(fluxZSum, expectations[interp - 1]);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
