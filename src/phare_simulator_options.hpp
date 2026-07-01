#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/def/phare_config.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE
{

namespace MHDOpts
{
    enum class TimeIntegratorType : uint8_t { Default, Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
    enum class ReconstructionType : uint8_t { Default, Constant, Linear, WENO3, WENOZ, MP5, count };
    enum class SlopeLimiterType : uint8_t { None, VanLeer, MinMod, count };
    enum class RiemannSolverType : uint8_t { Default, Rusanov, HLL, HLLD, count };
}; // namespace MHDOpts


struct SimOpts
{
    std::size_t dimension     = 1;
    std::size_t interp_order  = 1;
    std::size_t nbRefinedPart = core::defaultNbrRefinedParts(dimension, interp_order);

    core::LayoutMode layout_mode = core::LayoutMode::AoSMapped;
    AllocatorMode alloc_mode     = AllocatorMode::CPU;


    MHDOpts::TimeIntegratorType time_integrator_type = MHDOpts::TimeIntegratorType::Default;
    MHDOpts::ReconstructionType reconstruction_type  = MHDOpts::ReconstructionType::Default;
    MHDOpts::SlopeLimiterType slope_limiter_type     = MHDOpts::SlopeLimiterType::None;
    MHDOpts::RiemannSolverType riemann_solver_type   = MHDOpts::RiemannSolverType::Default;
    bool Hall                                        = false;
    bool Resistivity                                 = false;
    bool HyperResistivity                            = false;

    auto static constexpr make(std::size_t const dim, std::size_t const interp,
                               std::size_t const nbRefinedPart)
    {
        return SimOpts{.dimension = dim, .interp_order = interp, .nbRefinedPart = nbRefinedPart};
    }

    auto static constexpr make(std::size_t const dim, std::size_t const interp,
                               core::LayoutMode const layout_mode,
                               AllocatorMode const alloc_mode = AllocatorMode::CPU)
    {
        return SimOpts{dim, interp, core::defaultNbrRefinedParts(dim, interp), layout_mode,
                       alloc_mode};
    }

    auto static constexpr make(std::size_t const dim, std::size_t const interp,
                               std::size_t const nbRefinedPart, core::LayoutMode const layout_mode)
    {
        return SimOpts{dim, interp, nbRefinedPart, layout_mode, AllocatorMode::CPU};
    }

    auto static constexpr make(std::size_t const dim, std::size_t const interp,
                               std::size_t const nbRefinedPart, auto&&... args)
    {
        return SimOpts{
            dim, interp, nbRefinedPart, core::LayoutMode::AoSMapped, AllocatorMode::CPU, args...};
    }


    template<auto>
    struct Particles;
};

template<auto sopts>
struct SimOpts::Particles
{
    constexpr static core::ParticleArrayOptions opts{sopts.dimension, sopts.layout_mode,
                                                     core::StorageMode::VECTOR, sopts.alloc_mode};
    using value_type = core::ParticleArray<opts>;
};

} // namespace PHARE

#endif // PHARE_SIMULATOR_OPTIONS_HPP
