#ifndef PHARE_AMR_MESSENGERS_FUSER_OPS_STRATEGY_HPP
#define PHARE_AMR_MESSENGERS_FUSER_OPS_STRATEGY_HPP

// launches SAMRAI's staged kernel fusers around a schedule's copies: SAMRAI's schedules hand
//  fuseable transactions the fusers but never launch them, that is left to a ScheduleOpsStrategy.
//  FieldData/TensorFieldData::copyFuseable enqueue GPU_UNIFIED copies (mem::copy_fuseable), so
//  this is only attached to GPU_UNIFIED schedules, and launches only if something was enqueued

#include "core/def/phare_config.hpp"

#include <SAMRAI/tbox/ScheduleOpsStrategy.h>
#include <SAMRAI/tbox/StagedKernelFusers.h>

#include <cstddef>
#include <algorithm>

namespace PHARE::amr
{

class FuserOpsStrategy;

// SetEqual whose contiguous spans are enqueued in SAMRAI's stage 0 kernel fuser (GPU_UNIFIED with
//  RAJA) while a FuserOpsStrategy's schedule is copying/packing/unpacking, which launches them;
//  else (a schedule without the strategy) copied right away, see FieldData::copyFuseable
template<typename D>
struct SetEqualFuseable
{
    using value_type = D;

    void operator()(auto const& d0) { d = d0; }

    static void span(D* const dst, D const* const src, std::size_t const size)
    {
#if PHARE_HAVE_RAJA
        if (fusing())
        {
            // the staged enqueue: it flags StagedKernelFusers::enqueued(), which getFuser(0)->enqueue
            //  does not, so the strategy would never launch (they'd run at someone else's launch)
            SAMRAI::tbox::StagedKernelFusers::getInstance()->enqueue(
                0, 0, static_cast<int>(size),
                [=] SAMRAI_HOST_DEVICE(int const i) { dst[i] = src[i]; });
            return;
        }
#endif
        std::copy(src, src + size, dst);
    }

    static bool fusing();

    D& d;
};

class FuserOpsStrategy : public SAMRAI::tbox::ScheduleOpsStrategy
{
public:
    static FuserOpsStrategy* instance()
    {
        static FuserOpsStrategy strategy;
        return &strategy;
    }

    void preCommunicate() override {}
    void preCopy() override { fusing_ = true; }
    void prePack() override { fusing_ = true; }
    void preUnpack() override { fusing_ = true; }
    void postCommunicate() override {}
    void postCopy() override { launch_(); }
    void postPack() override { launch_(); }
    void postUnpack() override { launch_(); }
    bool deferMessageSend() override { return false; }

    // true between a pre and post hook of this strategy on this thread: only then are fuseable
    //  copies enqueued, as only then will they be launched
    static bool fusing() { return fusing_; }

private:
    static void launch_()
    {
        fusing_            = false;
        auto* const fusers = SAMRAI::tbox::StagedKernelFusers::getInstance();
        if (fusers->enqueued())
            fusers->launchAndCleanup();
    }

    static inline thread_local bool fusing_ = false;
};

template<typename D>
bool SetEqualFuseable<D>::fusing()
{
    return FuserOpsStrategy::fusing();
}

} // namespace PHARE::amr

#endif /* PHARE_AMR_MESSENGERS_FUSER_OPS_STRATEGY_HPP */
