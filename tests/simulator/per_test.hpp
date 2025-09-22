#ifndef PHARE_TEST_SIMULATOR_PER_TEST_HPP
#define PHARE_TEST_SIMULATOR_PER_TEST_HPP

#include "core/vector.hpp"
#include "phare/phare.hpp"
#include "core/def/phare_config.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "initializer/python_data_provider.hpp"

#include "simulator/simulator.hpp"

// #include "tests/core/data/field/test_field.hpp"

#include "gtest/gtest.h"


using SimOpts = PHARE::SimOpts;
using enum PHARE::core::LayoutMode;
using AllocatorMode = PHARE::AllocatorMode;


struct __attribute__((visibility("hidden"))) StaticIntepreter
{
    static std::shared_ptr<PHARE::initializer::PythonDataProvider> input;

    static StaticIntepreter& INSTANCE()
    {
        static StaticIntepreter i;
        return i;
    }
};
std::shared_ptr<PHARE::initializer::PythonDataProvider> StaticIntepreter::input{nullptr};


template<std::size_t dim>
struct HierarchyMaker
{
    HierarchyMaker(PHARE::initializer::PHAREDict& dict)
        : hierarchy{std::make_shared<PHARE::amr::DimHierarchy<dim>>(dict)}
    {
    }
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy;
};



template<auto opts>
struct SimulatorTestParam : private HierarchyMaker<opts.dimension>, public PHARE::Simulator<opts>
{
    static constexpr std::size_t dim = opts.dimension;


    using Simulator   = PHARE::Simulator<opts>;
    using PHARETypes  = PHARE::PHARE_Types<opts>;
    using Hierarchy   = PHARE::amr::Hierarchy;
    using HybridModel = PHARETypes::HybridModel_t;
    using MHDModel    = PHARETypes::MHDModel_t;
    using HierarchyMaker<dim>::hierarchy;

    auto& dict(std::string job_py)
    {
        auto& input = StaticIntepreter::INSTANCE().input;
        if (!input)
        {
            input = std::make_shared<PHARE::initializer::PythonDataProvider>(job_py);
            input->read();
        }
        SAMRAI::hier::VariableDatabase::getDatabase();
        return PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    }

    SimulatorTestParam(std::string job_py = "job")
        : HierarchyMaker<dim>{dict(job_py)}
        , Simulator{dict(job_py), this->hierarchy}
    {
        Simulator::initialize();
    }

    void reset_dman() { Simulator::reset_dman(); }
};

template<typename SimulatorParam>
struct SimulatorTest : public ::testing::Test
{
};

using Simulators = testing::Types<SimulatorTestParam<SimOpts::make(1, 1, 2)>>;
TYPED_TEST_SUITE(SimulatorTest, Simulators);


template<typename SimulatorParam>
struct Simulator1dTest : public ::testing::Test
{
};

// clang-format off
using Simulators1d = testing::Types<
    SimulatorTestParam<SimOpts::make(1, 1, 2)>, SimulatorTestParam<SimOpts::make(1, 1, 3)>,
    SimulatorTestParam<SimOpts::make(1, 2, 2)>, SimulatorTestParam<SimOpts::make(1, 2, 3)>,
    SimulatorTestParam<SimOpts::make(1, 2, 4)>, SimulatorTestParam<SimOpts::make(1, 3, 2)>,
    SimulatorTestParam<SimOpts::make(1, 3, 3)>, SimulatorTestParam<SimOpts::make(1, 3, 4)>,
    SimulatorTestParam<SimOpts::make(1, 3, 5)>

PHARE_WITH_MKN_GPU(
   ,SimulatorTestParam<SimOpts{1, 1, AoSTS, AllocatorMode::CPU}>
)
>;

TYPED_TEST_SUITE(Simulator1dTest, Simulators1d);


template<typename SimulatorParam>
struct Simulator2dTest : public ::testing::Test
{
};


using Simulators2d = testing::Types<
    SimulatorTestParam<SimOpts::make(2, 1, 4)>, SimulatorTestParam<SimOpts::make(2, 1, 5)>,
    SimulatorTestParam<SimOpts::make(2, 1, 8)>, SimulatorTestParam<SimOpts::make(2, 1, 9)>,
    SimulatorTestParam<SimOpts::make(2, 2, 4)>, SimulatorTestParam<SimOpts::make(2, 2, 5)>,
    SimulatorTestParam<SimOpts::make(2, 2, 8)>, SimulatorTestParam<SimOpts::make(2, 2, 9)>,
    SimulatorTestParam<SimOpts::make(2, 2, 16)>, SimulatorTestParam<SimOpts::make(2, 3, 4)>,
    SimulatorTestParam<SimOpts::make(2, 3, 5)>, SimulatorTestParam<SimOpts::make(2, 3, 8)>,
    SimulatorTestParam<SimOpts::make(2, 3, 25)>

PHARE_WITH_MKN_GPU(
   ,SimulatorTestParam<SimOpts{2, 1, AoSTS, AllocatorMode::CPU}>
)

>;
TYPED_TEST_SUITE(Simulator2dTest, Simulators2d);


template<typename Simulator>
struct Simulator3dTest : public ::testing::Test
{
};
using Simulator3d = testing::Types<
    SimulatorTestParam<SimOpts::make(3, 1, 6)>, SimulatorTestParam<SimOpts::make(3, 2, 6)>,
                       SimulatorTestParam<SimOpts::make(3, 3, 6)>
>;
TYPED_TEST_SUITE(Simulator3dTest, Simulator3d);

// clang-format on


#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */
