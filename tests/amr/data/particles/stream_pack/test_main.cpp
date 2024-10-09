
#include "core/def/phare_mpi.hpp"


#include <memory>
#include <cstdint>

#include "tests/amr/amr.hpp"

#include "core/logger.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/particles/particle_array_partitionner.hpp"

#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_data_factory.hpp"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "core/utilities/types.hpp"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;


template<std::size_t dim, typename ParticleArray>
struct AParticlesData
{
    static constexpr auto dimension = dim;
    static constexpr auto ghosts    = 1;
    using ParticleArray_t           = ParticleArray;
    using Particle_t                = typename ParticleArray_t::Particle_t;

    SAMRAI::tbox::Dimension amr_dimension{dim};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{amr_dimension, 0},
                                 SAMRAI::hier::Index{amr_dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{amr_dimension, 10},
                                   SAMRAI::hier::Index{amr_dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(amr_dimension)};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch destPatch{destDomain, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomain, patchDescriptor};

    ParticlesData<ParticleArray_t> destData{destDomain, ghost, "name"};
    ParticlesData<ParticleArray_t> sourceData{sourceDomain, ghost, "name"};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourceData.getGhostBox()};
    SAMRAI::hier::Box fillBox{destData.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(amr_dimension)};

    SAMRAI::hier::Transformation transformation{destDomain.lower() - sourceDomain.upper()
                                                - oneIndex};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};


    Particle_t particle;


    AParticlesData(std::array<int, dim> const& iCell)
    {
        particle.weight_ = 1.0;
        particle.charge_ = 1.0;
        particle.v_      = {1.0, 1.0, 1.0};
        particle.iCell() = iCell;
        sourceData.domainParticles.push_back(particle);
    }
};

template<std::size_t dim>
struct ACellMappedParticlesData : AParticlesData<dim, AoSMappedParticleArray<dim>>
{
    using Super = AParticlesData<dim, AoSMappedParticleArray<dim>>;

    ACellMappedParticlesData(std::array<int, dim> const& iCell)
        : Super{iCell}
    {
    }
};

template<std::size_t dim>
struct APartitionedParticlesData : AParticlesData<dim, AoSParticleArray<dim>>
{
    using Super    = AParticlesData<dim, AoSParticleArray<dim>>;
    using box_t    = PHARE::core::Box<int, dim>;
    using point_t  = Point<int, dim>;
    using iterator = typename std::vector<Particle<dim>>::iterator;
    using Super::cellOverlap;
    using Super::destDomain;
    using Super::sourceData;
    using Super::sourceDomain;
    using Super::transformation;


    auto static constexpr particle_ghosts = ghostWidthForParticles</*interp=*/1>();

    APartitionedParticlesData(std::array<int, dim> const& iCell)
        : Super{iCell}
    {
        auto destBox = destDomain;
        transformation.inverseTransform(destBox);
        src_neighbor_boxes.emplace_back(phare_box_from<dim>(destBox));

        PHARE_LOG_LINE_STR(phare_box_from<dim>(sourceDomain));
        PHARE_LOG_LINE_STR(src_neighbor_boxes.back());

        sourceData.set_iterators(partition<particle_ghosts>(
            sourceData.domainParticles, phare_box_from<dim>(sourceDomain), src_neighbor_boxes));
    }

    // std::vector<iterator> partition_iterators;
    std::vector<box_t> src_neighbor_boxes; //{phare_box_from<dim>(destDomain)};
    // std::vector<box_t> dst_neighbor_boxes{phare_box_from<dim>(sourceDomain)};
};

template<typename ParticlesData>
struct StreamPackTest : public ::testing::Test
{
};

using ParticlesDatas = testing::Types<ACellMappedParticlesData<
    1>                          /*, //
                                  // AParticlesData<2, AoSMappedParticleArray<2>>,
                                  // AParticlesData<3, AoSMappedParticleArray<3>>,
     APartitionedParticlesData<1> */ //
                                      >;
TYPED_TEST_SUITE(StreamPackTest, ParticlesDatas, );

TYPED_TEST(StreamPackTest, PreserveVelocityWhenPackStreamWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(15)};
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    // particle.iCell() = ConstArray<int, dim>(15);
    // sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].v_, Eq(particle.v_));
}




// TYPED_TEST(StreamPackTest, ShiftTheiCellWhenPackStreamWithPeriodics)
// {
//     using ParticlesData = TypeParam;
//     constexpr auto dim  = ParticlesData::dimension;

//     ParticlesData param{ConstArray<int, dim>(15)};
//     // auto& particle    = param.particle;
//     auto& sourceData  = param.sourceData;
//     auto& cellOverlap = param.cellOverlap;
//     auto& destData    = param.destData;

//     // particle.iCell() = ConstArray<int, dim>(15);
//     // sourceData.domainParticles.push_back(particle);

//     SAMRAI::tbox::MessageStream particlesWriteStream;

//     sourceData.packStream(particlesWriteStream, *cellOverlap);

//     SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
//                                                     SAMRAI::tbox::MessageStream::Read,
//                                                     particlesWriteStream.getBufferStart()};

//     destData.unpackStream(particlesReadStream, *cellOverlap);

//     // patch0 start at 0 , patch1 start at 10
//     // with periodics condition, we have 0 equivalent to 15
//     auto expectediCell = ConstArray<int, dim>(-1);
//     ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
//     ASSERT_THAT(destData.patchGhostParticles[0].iCell(), Eq(expectediCell));
// }



// TYPED_TEST(StreamPackTest, PackInTheCorrectBufferWithPeriodics)
// {
//     using ParticlesData = TypeParam;
//     constexpr auto dim  = ParticlesData::dimension;

//     ParticlesData param{ConstArray<int, dim>(15)};
//     // auto& particle    = param.particle;
//     auto& sourceData  = param.sourceData;
//     auto& cellOverlap = param.cellOverlap;
//     auto& destData    = param.destData;

//     // particle.iCell() = ;
//     // sourceData.domainParticles.push_back(particle);

//     SAMRAI::tbox::MessageStream particlesWriteStream;

//     sourceData.packStream(particlesWriteStream, *cellOverlap);

//     SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
//                                                     SAMRAI::tbox::MessageStream::Read,
//                                                     particlesWriteStream.getBufferStart()};

//     destData.unpackStream(particlesReadStream, *cellOverlap);

//     auto expectediCell = ConstArray<int, dim>(-1);

//     ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
//     ASSERT_THAT(destData.patchGhostParticles[0].iCell(), Eq(expectediCell));
// }



// TYPED_TEST(StreamPackTest,
//            PreserveParticleAttributesWhenPackingWithPeriodicsFromGhostSrcToDomainDest)
// {
//     using ParticlesData = TypeParam;
//     constexpr auto dim  = ParticlesData::dimension;

//     ParticlesData param{ConstArray<int, dim>(16)};
//     auto& particle    = param.particle;
//     auto& sourceData  = param.sourceData;
//     auto& cellOverlap = param.cellOverlap;
//     auto& destData    = param.destData;

//     // particle.iCell() = ;
//     // sourceData.domainParticles.push_back(particle);

//     SAMRAI::tbox::MessageStream particlesWriteStream;

//     sourceData.packStream(particlesWriteStream, *cellOverlap);

//     SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
//                                                     SAMRAI::tbox::MessageStream::Read,
//                                                     particlesWriteStream.getBufferStart()};

//     destData.unpackStream(particlesReadStream, *cellOverlap);

//     auto expectediCell = ConstArray<int, dim>(0);

// <<<<<<< HEAD
// auto expectediCell = ConstArray<int, dim>(0);

// EXPECT_THAT(destData.domainParticles[0].v, Eq(particle.v));
// EXPECT_THAT(destData.domainParticles[0].iCell, Eq(expectediCell));
// EXPECT_THAT(destData.domainParticles[0].delta, Eq(particle.delta));
// EXPECT_THAT(destData.domainParticles[0].weight, Eq(particle.weight));
// EXPECT_THAT(destData.domainParticles[0].charge, Eq(particle.charge));
// }
// =======
// //     EXPECT_THAT(destData.domainParticles[0].v_, Eq(particle.v_));
// //     EXPECT_THAT(destData.domainParticles[0].iCell_, Eq(expectediCell));
// //     EXPECT_THAT(destData.domainParticles[0].delta_, Eq(particle.delta_));
// //     EXPECT_THAT(destData.domainParticles[0].weight_, Eq(particle.weight_));
// //     EXPECT_THAT(destData.domainParticles[0].charge_, Eq(particle.charge_));
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].E_[0], particle.E_[0]);
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].E_[1], particle.E_[1]);
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].E_[2], particle.E_[2]);
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].B_[0], particle.B_[0]);
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].B_[1], particle.B_[1]);
// //     EXPECT_DOUBLE_EQ(destData.domainParticles[0].B_[2], particle.B_[2]);
// // }
// >>>>>>> 5e4f7fde (...)



// TYPED_TEST(StreamPackTest,
//            PreserveParticleAttributesWhenPackingWithPeriodicsFromDomainSrcToGhostDest)
// {
//     using ParticlesData = TypeParam;
//     constexpr auto dim  = ParticlesData::dimension;

//     ParticlesData param{ConstArray<int, dim>(15)};
//     auto& particle    = param.particle;
//     auto& sourceData  = param.sourceData;
//     auto& cellOverlap = param.cellOverlap;
//     auto& destData    = param.destData;

//     // particle.iCell() = ;
//     // sourceData.domainParticles.push_back(particle);

//     SAMRAI::tbox::MessageStream particlesWriteStream;

//     sourceData.packStream(particlesWriteStream, *cellOverlap);

//     SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
//                                                     SAMRAI::tbox::MessageStream::Read,
//                                                     particlesWriteStream.getBufferStart()};

//     destData.unpackStream(particlesReadStream, *cellOverlap);

//     auto expectediCell = ConstArray<int, dim>(-1);

// <<<<<<< HEAD
//     auto expectediCell = ConstArray<int, dim>(-1);

//     EXPECT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
//     EXPECT_THAT(destData.patchGhostParticles[0].iCell, Eq(expectediCell));
//     EXPECT_THAT(destData.patchGhostParticles[0].delta, Eq(particle.delta));
//     EXPECT_THAT(destData.patchGhostParticles[0].weight, Eq(particle.weight));
//     EXPECT_THAT(destData.patchGhostParticles[0].charge, Eq(particle.charge));
// }
// =======
// //     EXPECT_THAT(destData.patchGhostParticles[0].v_, Eq(particle.v_));
// //     EXPECT_THAT(destData.patchGhostParticles[0].iCell_, Eq(expectediCell));
// //     EXPECT_THAT(destData.patchGhostParticles[0].delta_, Eq(particle.delta_));
// //     EXPECT_THAT(destData.patchGhostParticles[0].weight_, Eq(particle.weight_));
// //     EXPECT_THAT(destData.patchGhostParticles[0].charge_, Eq(particle.charge_));
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].E_[0], particle.E_[0]);
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].E_[1], particle.E_[1]);
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].E_[2], particle.E_[2]);
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].B_[0], particle.B_[0]);
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].B_[1], particle.B_[1]);
// //     EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].B_[2], particle.B_[2]);
// // }
// >>>>>>> 5e4f7fde (...)



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::test::amr::SamraiLifeCycle life{argc, argv};
    return RUN_ALL_TESTS();
}
