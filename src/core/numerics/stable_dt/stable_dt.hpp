#ifndef PHARE_CORE_NUMERICS_STABLE_DT_STABLE_DT_HPP
#define PHARE_CORE_NUMERICS_STABLE_DT_STABLE_DT_HPP

#include "core/def.hpp"
#include "core/def/phare_config.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

#include "core/numerics/stable_dt/detail/detail.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <cstddef>

namespace PHARE::core
{
// CFL constraints as inverse time steps, maxed over a patch (or tile) domain
//  dt = cfl / max(inv_dt), the caller does the min over patches/tiles


// largest fn(cell, args...) over the domain box of field, fn must return a non negative value
template<typename Field, typename Fn, typename... Args>
double max_on_domain(auto const& layout, Field const& field, Fn const& fn, Args const&... args)
{
    auto const box = layout.domainBoxFor(field);

    if constexpr (Field::alloc_mode == AllocatorMode::CPU)
    {
        double max = 0;
        for_each_point(box, [&](auto const& bix) { max = std::max(max, fn(bix, args...)); });
        return max;
    }
    else
        return gpu::StableDt::max_on_box(fn, box, args...);
}


// pi * c_whistler, the hall whistler speed along a direction of mesh size dx
auto whistler_speed(auto const& dx, auto const& rho, auto const& BdotB) _PHARE_ALL_FN_
{
    return std::numbers::pi * compute_whistler_(1.0 / dx, rho, BdotB);
}


// B components projected onto the rho centering
template<typename GridLayout>
auto B_on_rho(auto const& Bx, auto const& By, auto const& Bz, auto const& index) _PHARE_ALL_FN_
{
    using Scalar = GridLayout::Scalar;
    return std::array{GridLayout::template project<Scalar::Bx, Scalar::rho>(Bx, index),
                      GridLayout::template project<Scalar::By, Scalar::rho>(By, index),
                      GridLayout::template project<Scalar::Bz, Scalar::rho>(Bz, index)};
}


// hybrid: whistler wave
template<typename GridLayout>
double max_inv_whistler_dt(GridLayout const& layout, auto const& n, auto const& B)
{
    auto constexpr static dim = GridLayout::dimension;

    auto const cell = [] _PHARE_ALL_FN_(auto const& bix, auto const& meshSize, auto const& n,
                                        auto const& Bx, auto const& By, auto const& Bz) {
        MeshIndex<dim> const index{bix};

        auto const rho          = n(index);
        auto const [bx, by, bz] = B_on_rho<GridLayout>(Bx, By, Bz, index);
        auto const BdotB        = bx * bx + by * by + bz * bz;

        double inv = 0;
        for (std::size_t d = 0; d < dim; ++d)
            inv += whistler_speed(meshSize[d], rho, BdotB) / meshSize[d];
        return inv;
    };

    return max_on_domain(layout, n, cell, layout.meshSize(), n, B(Component::X), B(Component::Y),
                         B(Component::Z));
}


// mhd: sum_d (|v_d| + c_fast_d + c_whistler_d) / dx_d, whistler only if hall
template<typename GridLayout>
double max_inv_mhd_wave_dt(GridLayout const& layout, double const gamma, bool const hall,
                           auto const& rho, auto const& rhoV, auto const& B, auto const& Etot)
{
    auto constexpr static dim = GridLayout::dimension;

    auto const cell = [] _PHARE_ALL_FN_(auto const& bix, auto const& meshSize, double const gamma,
                                        bool const hall, auto const& rho, auto const& rhoVx,
                                        auto const& rhoVy, auto const& rhoVz, auto const& Bx,
                                        auto const& By, auto const& Bz, auto const& Etot) {
        MeshIndex<dim> const index{bix};

        auto const r = rho(index);
        std::array const v{rhoVx(index) / r, rhoVy(index) / r, rhoVz(index) / r};
        auto const b     = B_on_rho<GridLayout>(Bx, By, Bz, index);
        auto const BdotB = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        auto const P     = eosEtotToP(gamma, r, v[0], v[1], v[2], b[0], b[1], b[2], Etot(index));

        double inv = 0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const cfast = compute_fast_magnetosonic_(gamma, r, b[d], BdotB, P);
            auto const cw    = hall ? whistler_speed(meshSize[d], r, BdotB) : 0.0;
            inv += (std::abs(v[d]) + cfast + cw) / meshSize[d];
        }
        return inv;
    };

    return max_on_domain(layout, rho, cell, layout.meshSize(), gamma, hall, rho, rhoV(Component::X),
                         rhoV(Component::Y), rhoV(Component::Z), B(Component::X), B(Component::Y),
                         B(Component::Z), Etot);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_NUMERICS_STABLE_DT_STABLE_DT_HPP */
