#ifndef PHARE_ION_UPDATER_DEF_HPP
#define PHARE_ION_UPDATER_DEF_HPP

#include "core/utilities/types.hpp"
#include <cstdint>

namespace PHARE::core
{
enum class UpdaterMode : std::uint16_t { domain_only = 0, all };
} // namespace PHARE::core


// ###############################################

#if PHARE_HAVE_MKN_GPU || defined(_MKN_WITH_MKN_KUL_)
#include "mkn/kul/os.hpp"

namespace PHARE::core::detail
{
auto static const timings_dir_str
    = get_env_as("PHARE_ASYNC_TIMES", std::string{".phare/async/multi_updater"});

static bool ion_updater_io_setup = []() {
    PHARE_LOG_LINE_SS("MAKE DIR!")
    mkn::kul::Dir timings{timings_dir_str};
    timings.mk();
    return true;
}();

} // namespace PHARE::core::detail

#endif // PHARE_HAVE_MKN_GPU

// ###############################################


#endif // ION_UPDATER_DEF_HPP
