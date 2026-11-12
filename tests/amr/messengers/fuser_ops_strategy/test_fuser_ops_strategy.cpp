// amr::SetEqualFuseable spans enqueued in SAMRAI's kernel fuser inside a FuserOpsStrategy
//  window are run (on the device for a CUDA SAMRAI) at its launch, outside one they copy at once

#include "core/def/phare_config.hpp"
#include "core/vector.hpp"
#include "amr/messengers/fuser_ops_strategy.hpp"

#include "tests/amr/amr.hpp"

#include "gtest/gtest.h"

#include <cstddef>

namespace PHARE::amr
{

using Vec_t = PHARE::Vector<double, AllocatorMode::GPU_UNIFIED>::vector_t;

struct Spans
{
    std::size_t static constexpr n_spans = 7, span = 1000;

    Spans()
        : src(n_spans * span)
        , dst(n_spans * span, -1.)
    {
        for (std::size_t i = 0; i < src.size(); ++i)
            src[i] = i;
    }

    void copy()
    {
        for (std::size_t s = 0; s < n_spans; ++s)
            SetEqualFuseable<double>::span(dst.data() + s * span, src.data() + s * span, span);
    }

    void expect_copied() const
    {
        for (std::size_t i = 0; i < src.size(); ++i)
            ASSERT_EQ(dst[i], src[i]) << i;
    }

    Vec_t src, dst;
};

TEST(FuserOpsStrategy, copies_outside_a_window_are_immediate)
{
    Spans spans;
    EXPECT_FALSE(FuserOpsStrategy::fusing());
    spans.copy();
    EXPECT_FALSE(SAMRAI::tbox::StagedKernelFusers::getInstance()->enqueued());
    spans.expect_copied();
}

TEST(FuserOpsStrategy, copies_in_a_window_run_at_its_launch)
{
    Spans spans;
    auto& strategy = *FuserOpsStrategy::instance();
    strategy.preCopy();
    EXPECT_TRUE(FuserOpsStrategy::fusing());
    spans.copy();
#if PHARE_HAVE_RAJA
    EXPECT_TRUE(SAMRAI::tbox::StagedKernelFusers::getInstance()->enqueued());
#endif
    strategy.postCopy(); // launch, synchronize, cleanup
    EXPECT_FALSE(FuserOpsStrategy::fusing());
    EXPECT_FALSE(SAMRAI::tbox::StagedKernelFusers::getInstance()->enqueued());
    spans.expect_copied();
}

} // namespace PHARE::amr


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::test::amr::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
