#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/def.hpp"

#include <cstdint>

namespace PHARE::core
{

enum class AmpereMode : std::uint16_t { ShrinkedGhost, GrownPhysical };

struct AmpereBox // structural type — usable as NTTP
{
    AmpereMode mode      = AmpereMode::ShrinkedGhost;
    std::uint32_t amount = 1;
};


// J = curl(B), one component at one point: curlc is that component of curl(layout, B)
//  where and how it is evaluated (boxes, tiles, GPU) is up to the caller,
//  see core/numerics/ampere/ampere_transformer.hpp
struct Ampere
{
    void operator()(auto const& ijk, auto& Jc, auto const& curlc) const _PHARE_ALL_FN_
    {
        Jc(ijk) = curlc({ijk});
    }
};

} // namespace PHARE::core

#endif
