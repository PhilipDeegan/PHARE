

#include "core/utilities/types.h"
#include "core/utilities/span.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/particle_initializers/particle_initializer_factory.h"
#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <vector>
#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::initializer;

using GridLayoutT    = GridLayout<GridLayoutImplYee<1, 1>>;
using ParticleArrayT = ParticleArray<1>;


namespace scalar_fn
{
double density(double /*x*/)
{
    return 1.;
}



double vx(double /*x*/)
{
    return 1.;
}


double vy(double /*x*/)
{
    return 1.;
}


double vz(double /*x*/)
{
    return 1.;
}


double vthx(double /*x*/)
{
    return 1.;
}


double vthy(double /*x*/)
{
    return 1.;
}



double vthz(double /*x*/)
{
    return 1.;
}

} // namespace scalar_fn

TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleInitializer)
{
    using namespace scalar_fn;
    PHARE::initializer::PHAREDict dict;
    dict["name"]            = std::string{"MaxwellianParticleInitializer"};
    dict["density"]         = static_cast<PHARE::initializer::ScalarFunction<1>>(density);
    dict["bulk_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);
    dict["bulk_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);
    dict["bulk_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);

    dict["thermal_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthx);
    dict["thermal_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthy);
    dict["thermal_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthz);

    dict["charge"]         = 1.;
    dict["nbrPartPerCell"] = int{100};
    dict["basis"]          = std::string{"Cartesian"};

    auto initializer = ParticleInitializerFactory<ParticleArrayT, GridLayoutT>::create(dict);
}

namespace vector_fn
{
template<typename T, typename SIZE = size_t>
class VectorSpan : public Span<T, SIZE>
{
    Span<T, SIZE> super(std::size_t size, T value)
    {
        vec = std::vector<T>(size, value);
        return Span<T, SIZE>{vec.data(), vec.size()};
    }

public:
    VectorSpan(std::size_t size, T value)
        : Span<T, SIZE>{super(size, value)}
    {
    }

private:
    std::vector<T> vec;
};

std::shared_ptr<PHARE::core::Span<double>> density(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}


std::shared_ptr<PHARE::core::Span<double>> vx(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

std::shared_ptr<PHARE::core::Span<double>> vy(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

std::shared_ptr<PHARE::core::Span<double>> vz(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

std::shared_ptr<PHARE::core::Span<double>> vthx(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

std::shared_ptr<PHARE::core::Span<double>> vthy(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}


std::shared_ptr<PHARE::core::Span<double>> vthz(std::vector<double> const& x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}
} // namespace vector_fn

TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleVectorInitializer)
{
    using namespace vector_fn;
    PHARE::initializer::PHAREDict dict;
    dict["name"]            = std::string{"MaxwellianParticleInitializer"};
    dict["fn_type"]         = std::string{"vector"};
    dict["density"]         = static_cast<PHARE::initializer::VectorFunction<1>>(density);
    dict["bulk_velocity_x"] = static_cast<PHARE::initializer::VectorFunction<1>>(vx);
    dict["bulk_velocity_y"] = static_cast<PHARE::initializer::VectorFunction<1>>(vx);
    dict["bulk_velocity_z"] = static_cast<PHARE::initializer::VectorFunction<1>>(vx);

    dict["thermal_velocity_x"] = static_cast<PHARE::initializer::VectorFunction<1>>(vthx);
    dict["thermal_velocity_y"] = static_cast<PHARE::initializer::VectorFunction<1>>(vthy);
    dict["thermal_velocity_z"] = static_cast<PHARE::initializer::VectorFunction<1>>(vthz);

    dict["charge"]         = 1.;
    dict["nbrPartPerCell"] = int{100};
    dict["basis"]          = std::string{"Cartesian"};

    auto initializer = ParticleInitializerFactory<ParticleArrayT, GridLayoutT>::create(dict);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
