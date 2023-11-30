#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/index/index.hpp"

#include "phare_core.hpp"


using namespace PHARE::core;


template<int dim, int interp>
class NDlayout
{
    NDlayout() {}

    using nDL = GridLayout<GridLayoutImplYee<dim, interp>>;

public:
    static nDL create()
    {
        if constexpr (dim == 1)
        {
            return {{{0.1}}, {{50}}, {0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{0.1, 0.2}}, {{50, 40}}, {0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{0.1, 0.2, 0.3}}, {{50, 40, 30}}, {0., 0., 0.}};
        }
    }
};




PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["resistivity"]       = 1.0;
    dict["hyper_resistivity"] = 0.01;

    return dict;
}



template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct OhmTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee  = GridLayout<GridLayoutImplYee<dim, interp>>;
    using Grid_t   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    GridYee layout = NDlayout<dim, interp>::create();

    Grid_t n;
    Grid_t Vx;
    Grid_t Vy;
    Grid_t Vz;
    Grid_t P;
    Grid_t Bx;
    Grid_t By;
    Grid_t Bz;
    Grid_t Jx;
    Grid_t Jy;
    Grid_t Jz;
    Grid_t Exnew;
    Grid_t Eynew;
    Grid_t Eznew;
    VecField<Grid_t, HybridQuantity> V;
    VecField<Grid_t, HybridQuantity> B;
    VecField<Grid_t, HybridQuantity> J;
    VecField<Grid_t, HybridQuantity> Enew;
    Ohm<GridYee> ohm;

    OhmTest()
        : n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
        , Vx{"Vx", HybridQuantity::Scalar::Vx, layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"Vy", HybridQuantity::Scalar::Vy, layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"Vz", HybridQuantity::Scalar::Vz, layout.allocSize(HybridQuantity::Scalar::Vz)}
        , P{"P", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Eynew{"Eynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Eznew{"Eznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , V{"V", HybridQuantity::Vector::V}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
        , Enew{"Enew", HybridQuantity::Vector::E}
        , ohm{createDict()}
    {
        V.setBuffer("V_x", &Vx);
        V.setBuffer("V_y", &Vy);
        V.setBuffer("V_z", &Vz);
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
        Enew.setBuffer("Enew_x", &Exnew);
        Enew.setBuffer("Enew_y", &Eynew);
        Enew.setBuffer("Enew_z", &Eznew);

        if constexpr (dim == 1)
        {
            auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
            auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
            auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
            auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

            for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
            {
                auto point = this->layout.fieldNodeCoordinates(n, Point<double, 1>{0.}, ix);

                n(ix)  = std::cosh(0.5 * point[0]);
                Vx(ix) = std::sinh(0.2 * point[0]);
                Vy(ix) = std::sinh(0.3 * point[0]);
                Vz(ix) = std::sinh(0.4 * point[0]);
                P(ix)  = std::cosh(0.5 * point[0]);
                Bx(ix) = std::cosh(0.2 * point[0]);
                Jy(ix) = std::tanh(0.3 * point[0]);
                Jz(ix) = std::tanh(0.4 * point[0]);
            }

            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                auto point = this->layout.fieldNodeCoordinates(Bz, Point<double, 1>{0.}, ix);

                By(ix) = std::cosh(0.3 * point[0]);
                Bz(ix) = std::cosh(0.4 * point[0]);
                Jx(ix) = std::tanh(0.2 * point[0]);
            }
        }

        if constexpr (dim == 2)
        {
            auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
            auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
            auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
            auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
            auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
            auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
            auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
            auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

            for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
            {
                for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
                {
                    auto point
                        = this->layout.fieldNodeCoordinates(n, Point<double, 2>{0., 0.}, ix, iy);

                    n(ix, iy)  = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1]);
                    Vx(ix, iy) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1]);
                    Vy(ix, iy) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1]);
                    Vz(ix, iy) = std::sinh(0.4 * point[0]) * std::sinh(0.4 * point[1]);
                    P(ix, iy)  = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1]);
                    Jz(ix, iy) = std::tanh(0.4 * point[0]) * std::tanh(0.4 * point[1]);
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    auto point
                        = this->layout.fieldNodeCoordinates(Bx, Point<double, 2>{0., 0.}, ix, iy);

                    Bx(ix, iy) = std::cosh(0.2 * point[0]) * std::cosh(0.2 * point[1]);
                    Jy(ix, iy) = std::tanh(0.3 * point[0]) * std::tanh(0.3 * point[1]);
                }
            }
            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
                {
                    auto point
                        = this->layout.fieldNodeCoordinates(Jx, Point<double, 2>{0., 0.}, ix, iy);

                    By(ix, iy) = std::cosh(0.3 * point[0]) * std::cosh(0.3 * point[1]);
                    Jx(ix, iy) = std::tanh(0.2 * point[0]) * std::tanh(0.2 * point[1]);
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    auto point
                        = this->layout.fieldNodeCoordinates(Bz, Point<double, 2>{0., 0.}, ix, iy);

                    Bz(ix, iy) = std::cosh(0.4 * point[0]) * std::cosh(0.4 * point[1]);
                }
            }
        }

        if constexpr (dim == 3)
        {
            auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
            auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
            auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
            auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
            auto gsi_p_Z = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Z);
            auto gei_p_Z = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Z);
            auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
            auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
            auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
            auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);
            auto gsi_d_Z = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Z);
            auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

            for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
            {
                for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            n, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        n(ix, iy, iz) = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1])
                                        * std::cosh(0.5 * point[2]);
                        Vx(ix, iy, iz) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1])
                                         * std::sinh(0.2 * point[2]);
                        Vy(ix, iy, iz) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1])
                                         * std::sinh(0.3 * point[2]);
                        Vz(ix, iy, iz) = std::sinh(0.4 * point[0]) * std::sinh(0.4 * point[1])
                                         * std::sinh(0.4 * point[2]);
                        P(ix, iy, iz) = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1])
                                        * std::cosh(0.5 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Jz, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        Jz(ix, iy, iz) = std::tanh(0.4 * point[0]) * std::tanh(0.4 * point[1])
                                         * std::tanh(0.4 * point[2]);
                    }
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Jy, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        Jy(ix, iy, iz) = std::tanh(0.3 * point[0]) * std::tanh(0.3 * point[1])
                                         * std::tanh(0.3 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Bx, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        Bx(ix, iy, iz) = std::cosh(0.2 * point[0]) * std::cosh(0.2 * point[1])
                                         * std::cosh(0.2 * point[2]);
                    }
                }
            }
            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Jx, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        Jx(ix, iy, iz) = std::tanh(0.2 * point[0]) * std::tanh(0.2 * point[1])
                                         * std::tanh(0.2 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            By, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        By(ix, iy, iz) = std::cosh(0.3 * point[0]) * std::cosh(0.3 * point[1])
                                         * std::cosh(0.3 * point[2]);
                    }
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Bz, Point<double, 3>{0., 0., 0.}, ix, iy, iz);

                        Bz(ix, iy, iz) = std::cosh(0.4 * point[0]) * std::cosh(0.4 * point[1])
                                         * std::cosh(0.4 * point[2]);
                    }
                }
            }
        }
    }

    ~OhmTest() {}
};



using OhmTupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;

TYPED_TEST_SUITE(OhmTest, OhmTupleInfos);




TYPED_TEST(OhmTest, ThatOhmHasCtorWithDict)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    Ohm<GridYee> ohm(createDict());
}




TYPED_TEST(OhmTest, ShouldBeGivenAGridLayoutPointerToBeOperational)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto layout = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    // this->ohm.setLayout(layout.get());
    EXPECT_ANY_THROW(
        this->ohm(this->n, this->V, this->P, this->B, this->J,
                  this->Enew)); // because the grid layout is not set (TODO issue #3392)
}


std::vector<double> read(std::string filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}



TYPED_TEST(OhmTest, ThatElectricFieldIsOkFromOhmsLaw)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    std::string prefix{"ohm"};

    auto filenameX
        = prefix + "x_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto filenameY
        = prefix + "y_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto filenameZ
        = prefix + "z_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto expected_ohmX = read(filenameX);
    auto expected_ohmY = read(filenameY);
    auto expected_ohmZ = read(filenameZ);

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    auto layout   = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    this->ohm.setLayout(layout.get());
    this->ohm(this->n, this->V, this->P, this->B, this->J, this->Enew);

    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(this->Exnew, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(this->Exnew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(this->Exnew(ix), ::testing::DoubleNear((expected_ohmX[ix]), 1e-12));
        }

        psi_X = this->layout.physicalStartIndex(this->Eynew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eynew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(this->Eynew(ix), ::testing::DoubleNear((expected_ohmY[ix]), 1e-12));
        }

        psi_X = this->layout.physicalStartIndex(this->Eznew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eznew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(this->Eznew(ix), ::testing::DoubleNear((expected_ohmZ[ix]), 1e-12));
        }
    }

    if constexpr (dim == 2)
    {
        auto psi_X = this->layout.physicalStartIndex(this->Exnew, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(this->Exnew, Direction::X);
        auto psi_Y = this->layout.physicalStartIndex(this->Exnew, Direction::Y);
        auto pei_Y = this->layout.physicalEndIndex(this->Exnew, Direction::Y);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ex);
                auto index_ = ix * nPts_[1] + iy;
                EXPECT_THAT(this->Exnew(ix, iy),
                            ::testing::DoubleNear((expected_ohmX[index_]), 1e-12));
            }
        }

        psi_X = this->layout.physicalStartIndex(this->Eynew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eynew, Direction::X);
        psi_Y = this->layout.physicalStartIndex(this->Eynew, Direction::Y);
        pei_Y = this->layout.physicalEndIndex(this->Eynew, Direction::Y);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ey);
                auto index_ = ix * nPts_[1] + iy;
                EXPECT_THAT(this->Eynew(ix, iy),
                            ::testing::DoubleNear((expected_ohmY[index_]), 1e-12));
            }
        }

        psi_X = this->layout.physicalStartIndex(this->Eznew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eznew, Direction::X);
        psi_Y = this->layout.physicalStartIndex(this->Eznew, Direction::Y);
        pei_Y = this->layout.physicalEndIndex(this->Eznew, Direction::Y);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ez);
                auto index_ = ix * nPts_[1] + iy;
                EXPECT_THAT(this->Eznew(ix, iy),
                            ::testing::DoubleNear((expected_ohmZ[index_]), 1e-12));
            }
        }
    }

    if constexpr (dim == 3)
    {
        auto psi_X = this->layout.physicalStartIndex(this->Exnew, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(this->Exnew, Direction::X);
        auto psi_Y = this->layout.physicalStartIndex(this->Exnew, Direction::Y);
        auto pei_Y = this->layout.physicalEndIndex(this->Exnew, Direction::Y);
        auto psi_Z = this->layout.physicalStartIndex(this->Exnew, Direction::Z);
        auto pei_Z = this->layout.physicalEndIndex(this->Exnew, Direction::Z);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                {
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ex);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(this->Exnew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmX[index_]), 1e-10));
                }
            }
        }

        psi_X = this->layout.physicalStartIndex(this->Eynew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eynew, Direction::X);
        psi_Y = this->layout.physicalStartIndex(this->Eynew, Direction::Y);
        pei_Y = this->layout.physicalEndIndex(this->Eynew, Direction::Y);
        psi_Z = this->layout.physicalStartIndex(this->Eynew, Direction::Z);
        pei_Z = this->layout.physicalEndIndex(this->Eynew, Direction::Z);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                {
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ey);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(this->Eynew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmY[index_]), 1e-10));
                }
            }
        }

        psi_X = this->layout.physicalStartIndex(this->Eznew, Direction::X);
        pei_X = this->layout.physicalEndIndex(this->Eznew, Direction::X);
        psi_Y = this->layout.physicalStartIndex(this->Eznew, Direction::Y);
        pei_Y = this->layout.physicalEndIndex(this->Eznew, Direction::Y);
        psi_Z = this->layout.physicalStartIndex(this->Eznew, Direction::Z);
        pei_Z = this->layout.physicalEndIndex(this->Eznew, Direction::Z);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                {
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ez);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(this->Eznew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmZ[index_]), 1e-10));
                }
            }
        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
