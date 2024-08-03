
#include "phare_core.hpp"

//#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/particles_data.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "tests/core/utilities/box/test_box_fixtures.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;
using namespace PHARE::amr;


std::size_t constexpr static dim = 3; // only 3d
using Box_t                      = PHARE::core::Box<int, dim>;
using FileParticles              = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoS>>;

auto file_name_for_box(Box_t const& box)
{
    std::stringstream ss;
    ss << "particles." << std::hex << std::hash<std::string>()(std::to_string(box)) << ".bin";
    return ss.str();
}

struct StaticParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = 1;
    std::size_t constexpr static ppc         = 100;
    std::size_t constexpr static box_size    = 50;

    static StaticParticlesDataFixture instance()
    {
        static StaticParticlesDataFixture i;
        return i;
    }

    StaticParticlesDataFixture()
        : box_tuple{get_some_boxes(box_size)}
    {
        auto write_particles = [](auto const& box) {
            FileParticles particles;
            add_particles_in(particles, box, ppc);
            disperse(particles, 1333337);
            write_raw_to_file(particles, file_name_for_box(box));
        };
        write_particles(middle());
    }

    auto b() const { return true; }
    auto& middle() const { return std::get<0>(box_tuple); }
    auto& neighbours() const { return std::get<1>(box_tuple); }
    std::tuple<Box_t, std::vector<Box_t>> const box_tuple;
};

static auto is_initialized = StaticParticlesDataFixture::instance().b();

// PP = per particle
template<typename Particles_>
struct ParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = StaticParticlesDataFixture::ghost_cells;
    std::size_t constexpr static ppc         = StaticParticlesDataFixture::ppc;
    std::size_t constexpr static box_size    = StaticParticlesDataFixture::box_size;

    SAMRAI::tbox::Dimension static inline const dimension{dim};
    SAMRAI::hier::BlockId static inline const blockId{0};
    SAMRAI::hier::IntVector static inline const ghost{dimension, ghost_cells};

    using Particles_t     = Particles_;
    using ParticlesData_t = ParticlesData<Particles_t>;

    ParticlesDataFixture()
        : middle{samrai_box_from(StaticParticlesDataFixture::instance().middle()), ghost}
        , neighbours{generate_from(
              [&](auto& box) {
                  return std::make_unique<ParticlesData_t>(samrai_box_from(box), ghost);
              },
              StaticParticlesDataFixture::instance().neighbours())}
    {
        for (auto const& particle : read_raw_from_file<FileParticles>(
                 file_name_for_box(phare_box_from<dim>(middle.getBox()))))
            middle.domainParticles.push_back(particle);
    }

    // auto static overlap(ParticlesDatas_t const& neighbour)
    // {
    //     std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
    //         std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};
    //     std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
    //         std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};

    //     SAMRAI::hier::Box srcMask{sourcePdat.getGhostBox()};
    //     SAMRAI::hier::Box fillBox{destPdat.getGhostBox()};

    //     bool overwriteInterior{true};
    //     SAMRAI::hier::Transformation transformation{destPdat.getGhostBox().lower()
    //                                                 - sourceDomain.upper()};

    //     return std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
    //         std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
    //             *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};
    // }

    ParticlesData_t middle;
    std::vector<std::unique_ptr<ParticlesData_t>> neighbours;
};


template<typename ParticlesData>
struct ParticlesDataTest : public ::testing::Test, public ParticlesData
{
};


using ParticlesDatas = testing::Types< //
    ParticlesDataFixture<ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoS>>>,
    ParticlesDataFixture<ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoSMapped>>>,
    ParticlesDataFixture<ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::SoA>>>,
    ParticlesDataFixture<
        ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoS, StorageMode::VECTOR,
                                                  PHARE::AllocatorMode::GPU_UNIFIED>>>,
    ParticlesDataFixture<
        ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::SoA, StorageMode::VECTOR,
                                                  PHARE::AllocatorMode::GPU_UNIFIED>>>

    >;

TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas);

TYPED_TEST(ParticlesDataTest, copyWorks)
{
    using Test = TypeParam;
    auto start = now_in_microseconds();
    for (auto& pdata : this->neighbours)
        this->middle.copy(*pdata);
    auto end = now_in_microseconds();
    PHARE_LOG_LINE_STR("time " << end - start);
    std::size_t expected
        = (std::pow(Test::box_size + (Test::ghost_cells * 2), dim) - std::pow(Test::box_size, dim))
          * Test::ppc;
    EXPECT_EQ(expected, this->middle.patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, copyWithTransformWorks)
{
    // using Test = TypeParam;
    // auto start = now_in_microseconds();
    // for (auto& pdata : this->neighbours)
    //     this->middle.copy(*pdata);
    // auto end = now_in_microseconds();
    // PHARE_LOG_LINE_STR("time " << end - start);
    // std::size_t expected
    //     = (std::pow(Test::box_size + (Test::ghost_cells * 2), dim) - std::pow(Test::box_size,
    //     dim))
    //       * Test::ppc;
    // EXPECT_EQ(expected, this->middle.patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packWorks)
{
    // using Test = TypeParam;
    // auto start = now_in_microseconds();
    // for (auto& pdata : this->neighbours)
    //     this->middle.copy(*pdata);
    // auto end = now_in_microseconds();
    // PHARE_LOG_LINE_STR("time " << end - start);
    // std::size_t expected
    //     = (std::pow(Test::box_size + (Test::ghost_cells * 2), dim) - std::pow(Test::box_size,
    //     dim))
    //       * Test::ppc;
    // EXPECT_EQ(expected, this->middle.patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packWithTransformWorks)
{
    // using Test = TypeParam;
    // auto start = now_in_microseconds();
    // for (auto& pdata : this->neighbours)
    //     this->middle.copy(*pdata);
    // auto end = now_in_microseconds();
    // PHARE_LOG_LINE_STR("time " << end - start);
    // std::size_t expected
    //     = (std::pow(Test::box_size + (Test::ghost_cells * 2), dim) - std::pow(Test::box_size,
    //     dim))
    //       * Test::ppc;
    // EXPECT_EQ(expected, this->middle.patchGhostParticles.size());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
