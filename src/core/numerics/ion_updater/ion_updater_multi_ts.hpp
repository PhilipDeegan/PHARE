#ifndef PHARE_ION_UPDATER_MULTI_TS_HPP
#define PHARE_ION_UPDATER_MULTI_TS_HPP

#include "ion_updater_def.hpp"

#include "core/numerics/pusher/multi_boris.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

#if PHARE_HAVE_MKN_GPU
#include "mkn/kul/os.hpp"
#endif // PHARE_HAVE_MKN_GPU


namespace PHARE::core::detail
{
// auto static const timings_dir_str
//     = get_env_as("PHARE_ASYNC_TIMES", std::string{".phare/async/multi_updater"});

static bool ion_updater_ts_setup = []() {
    mkn::kul::Dir timings{timings_dir_str};
    timings.mk();
    return true;
}();

} // namespace PHARE::core::detail

#include "core/numerics/ion_updater/ts_impls/entry.hpp"




#endif // ION_UPDATER_HPP
