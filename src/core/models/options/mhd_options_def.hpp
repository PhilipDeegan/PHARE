#ifndef PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP
#define PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP

#include <cstdint>

namespace PHARE
{
namespace MHDOpts
{
    enum class TimeIntegratorType : uint8_t { Default, Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
    enum class ReconstructionType : uint8_t { Default, Constant, Linear, WENO3, WENOZ, MP5, count };
    enum class SlopeLimiterType : uint8_t { None, VanLeer, MinMod, count };
    enum class RiemannSolverType : uint8_t { Default, Rusanov, HLL, HLLD, count };
}; // namespace MHDOpts

} // namespace PHARE

#endif // PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP
