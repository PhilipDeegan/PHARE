#ifndef PHARE_CORE_PUSHER_BORIS_BASICS_HPP
#define PHARE_CORE_PUSHER_BORIS_BASICS_HPP

#include "core/utilities/types.hpp"

#include <tuple>
#include <iterator>
#include <utility>

namespace PHARE::core::boris
{



template<auto alloc, typename Particle_t, typename Float, std::size_t dim>
auto static advance(Particle_t& p, std::array<Float, dim> const& halfDtOverDl) _PHARE_ALL_FN_
{
    auto newCell = p.iCell();

    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        Float delta = p.delta()[iDim] + static_cast<Float>(halfDtOverDl[iDim] * p.v()[iDim]);
        if constexpr (alloc == AllocatorMode::CPU)
            if (std::abs(delta) > 2)
                throw_runtime_error("Error, particle moves more than 1 cell, delta >2");

        Float iCell     = std::floor(delta);
        p.delta()[iDim] = delta - iCell;
        newCell[iDim]   = static_cast<int>(iCell + p.iCell()[iDim]);
    }

    return newCell;
}


template<typename Particle, typename ParticleEB, typename Float>
void accelerate(Particle& p, ParticleEB const& eb, Float const& dto2m) _PHARE_ALL_FN_
{
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    Float const coef1 = p.charge() * dto2m;

    // We now apply the 3 steps of the BORIS PUSHER
    // 1st half push of the electric field
    Float const velx1 = p.v()[0] + coef1 * pEx;
    Float const vely1 = p.v()[1] + coef1 * pEy;
    Float const velz1 = p.v()[2] + coef1 * pEz;

    // preparing variables for magnetic rotation
    Float const rx = coef1 * pBx;
    Float const ry = coef1 * pBy;
    Float const rz = coef1 * pBz;

    Float const rx2  = rx * rx;
    Float const ry2  = ry * ry;
    Float const rz2  = rz * rz;
    Float const rxry = rx * ry;
    Float const rxrz = rx * rz;
    Float const ryrz = ry * rz;

    Float const invDet = one / (one + rx2 + ry2 + rz2);

    // preparing rotation matrix due to the magnetic field
    // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
    Float const mxx = one + rx2 - ry2 - rz2;
    Float const mxy = two * (rxry + rz);
    Float const mxz = two * (rxrz - ry);

    Float const myx = two * (rxry - rz);
    Float const myy = one + ry2 - rx2 - rz2;
    Float const myz = two * (ryrz + rx);

    Float const mzx = two * (rxrz + ry);
    Float const mzy = two * (ryrz - rx);
    Float const mzz = one + rz2 - rx2 - ry2;

    // magnetic rotation
    Float const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
    Float const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
    Float const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;

    // 2nd half push of the electric field / Update particle velocity
    p.v()[0] = velx2 + coef1 * pEx;
    p.v()[1] = vely2 + coef1 * pEy;
    p.v()[2] = velz2 + coef1 * pEz;
}


template<typename Particles_t, typename... Args>
auto advance_noavx(Particles_t& ps, Args const&... args) _PHARE_ALL_FN_
{
    auto const& [pidx, halfDtOverDl] = std::forward_as_tuple(args...);
    using Float                      = std::decay_t<decltype(halfDtOverDl)>::value_type;

    auto newCell = ps.iCell(pidx);

    for (std::size_t iDim = 0; iDim < Particles_t::dimension; ++iDim)
    {
        auto& deldim     = ps.delta()[iDim][pidx];
        auto const delta = deldim + static_cast<Float>(halfDtOverDl[iDim] * ps.v()[iDim][pidx]);
        if constexpr (Particles_t::alloc_mode == AllocatorMode::CPU)
            if (std::abs(delta) > 2)
                throw_runtime_error("Error, particle moves more than 1 cell, delta >2");

        auto iCell = std::floor(delta);
        deldim     = delta - iCell;
        newCell[iDim] += static_cast<int>(iCell);
    }

    ps.iCell(pidx) = newCell;
}

template<typename Particles_t, typename... Args>
void accelerate_noavx(Particles_t& ps, Args&... args) _PHARE_ALL_FN_
{
    auto const& [pidx, em, interp, layout, halfDtOverDl, dto2m] = std::forward_as_tuple(args...);
    using Float                = std::decay_t<decltype(halfDtOverDl)>::value_type;
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto const& eb        = interp.m2p(ps.begin() + pidx, em, layout);
    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    Float const coef1 = ps.charge()[pidx] * dto2m;

    // We now apply the 3 steps of the BORIS PUSHER
    // 1st half push of the electric field
    Float const velx1 = ps.v()[0][pidx] + coef1 * pEx;
    Float const vely1 = ps.v()[1][pidx] + coef1 * pEy;
    Float const velz1 = ps.v()[2][pidx] + coef1 * pEz;

    // preparing variables for magnetic rotation
    Float const rx = coef1 * pBx;
    Float const ry = coef1 * pBy;
    Float const rz = coef1 * pBz;

    Float const rx2  = rx * rx;
    Float const ry2  = ry * ry;
    Float const rz2  = rz * rz;
    Float const rxry = rx * ry;
    Float const rxrz = rx * rz;
    Float const ryrz = ry * rz;

    Float const invDet = one / (one + rx2 + ry2 + rz2);

    // preparing rotation matrix due to the magnetic field
    // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
    Float const mxx = one + rx2 - ry2 - rz2;
    Float const mxy = two * (rxry + rz);
    Float const mxz = two * (rxrz - ry);

    Float const myx = two * (rxry - rz);
    Float const myy = one + ry2 - rx2 - rz2;
    Float const myz = two * (ryrz + rx);

    Float const mzx = two * (rxrz + ry);
    Float const mzy = two * (ryrz - rx);
    Float const mzz = one + rz2 - rx2 - ry2;

    // magnetic rotation
    Float const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
    Float const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
    Float const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;

    // 2nd half push of the electric field / Update particle velocity
    ps.v()[0][pidx] = velx2 + coef1 * pEx;
    ps.v()[1][pidx] = vely2 + coef1 * pEy;
    ps.v()[2][pidx] = velz2 + coef1 * pEz;
}


#if PHARE_HAVE_MKN_AVX



template<typename Particles_t, typename... Args>
auto advance_avx(Particles_t& ps, Args const&... args) _PHARE_ALL_FN_
{
    auto const& [pidx, halfDtOverDl] = std::forward_as_tuple(args...);
    using Float                      = std::decay_t<decltype(halfDtOverDl)>::value_type;

    constexpr auto N = mkn::avx::Span<Float>::N;

    for (std::size_t iDim = 0; iDim < Particles_t::dimension; ++iDim)
    {
        auto delta = mkn::avx::make_span(ps.delta()[iDim], pidx, N);
        {
            mkn::avx::Array<Float, N> newdelta{halfDtOverDl[iDim]};
            newdelta *= mkn::avx::make_span(ps.v()[iDim], pidx, N);
            newdelta += delta;
            delta = newdelta;
        }

        for (std::size_t i = 0; i < N; ++i)
        {
            auto iCell = std::floor(delta[i]);
            delta[i] -= iCell;
            ps.iCell(pidx + i)[iDim] += static_cast<int>(iCell);
        }
    }
}


template<typename Particles_t, typename... Args>
void accelerate_avx(Particles_t& ps, Args&... args) _PHARE_ALL_FN_
{
    auto const& [pidx, em, interp, layout, halfDtOverDl, dto2m] = std::forward_as_tuple(args...);
    using Float = std::decay_t<decltype(halfDtOverDl)>::value_type;

    constexpr auto N = mkn::avx::Span<Float>::N;
    static_assert(N == 8);
    static_assert(mkn::avx::Options::ALIGN() == 64);

    using Arr = mkn::avx::Array<Float, N>;

    Arr const one{1};
    Arr const two{2};

    auto const ebs = interp.template m2p_avx<N, mkn::avx::Array>(layout, ps, em, pidx);
    auto const& [pEx, pEy, pEz, pBx, pBy, pBz] = ebs;

    auto v0 = mkn::avx::make_span<N>(ps.v_[0], pidx);
    auto v1 = mkn::avx::make_span<N>(ps.v_[1], pidx);
    auto v2 = mkn::avx::make_span<N>(ps.v_[2], pidx);

    auto const charge = mkn::avx::make_span<N>(ps.charge_, pidx);
    auto const coef1  = Arr::FROM([&](auto const e) { return e * dto2m; }, charge);

    auto const vel = [&](auto const& v, auto const& p) {
        auto vel = coef1;
        vel *= p;
        vel += v;
        return vel;
    };

    auto const velx1 = vel(v0, pEx);
    auto const vely1 = vel(v1, pEy);
    auto const velz1 = vel(v2, pEz);


    auto const rx = coef1 * pBx;
    auto const ry = coef1 * pBy;
    auto const rz = coef1 * pBz;

    auto const rx2  = rx * rx;
    auto const ry2  = ry * ry;
    auto const rz2  = rz * rz;
    auto const rxry = rx * ry;
    auto const rxrz = rx * rz;
    auto const ryrz = ry * rz;

    // auto const invDet = one / (one + rx2 + ry2 + rz2);
    auto const _invDet = [&]() {
        auto det = one;
        det += rx2;
        det += ry2;
        det += rz2;
        return one / det;
    };
    auto const invDet = _invDet();


    // preparing rotation matrix due to the magnetic field

    auto const msubr2 = [&](auto const&... args) {
        auto const& [r20, r21, r22] = std::forward_as_tuple(args...);
        auto m                      = one;
        m += r20;
        m -= r21;
        m -= r22;
        return m;
    };
    auto const mmuladdr2 = [&](auto const& rr, auto const& r) {
        auto m = rr;
        m += r;
        m *= two;
        return m;
    };
    auto const mmulsubr2 = [&](auto const& rr, auto const& r) {
        auto m = rr;
        m -= r;
        m *= two;
        return m;
    };

    auto const mxx = msubr2(rx2, ry2, rz2);
    auto const mxy = mmuladdr2(rxry, rz);
    auto const mxz = mmulsubr2(rxrz, ry);

    auto const myx = mmulsubr2(rxry, rz);
    auto const myy = msubr2(ry2, rx2, rz2);
    auto const myz = mmuladdr2(ryrz, rx);

    auto const mzx = mmuladdr2(rxrz, ry);
    auto const mzy = mmulsubr2(ryrz, rx);
    auto const mzz = msubr2(rz2, rx2, ry2);


    // magnetic rotation
    auto const vel2 = [&](auto const&... args) {
        auto const& [x, y, z] = std::forward_as_tuple(args...);
        mkn::avx::Array<Float, N> vel{0};

        auto velx = velx1;
        auto vely = vely1;
        auto velz = velz1;

        velx *= x;
        vely *= y;
        velz *= z;

        vel += velx;
        vel += vely;
        vel += velz;

        return vel * invDet;
    };

    auto const velx2 = vel2(mxx, mxy, mxz);
    auto const vely2 = vel2(myx, myy, myz);
    auto const velz2 = vel2(mzx, mzy, mzz);

    // 2nd half push of the electric field / Update particle velocity
    auto const setvel = [&](auto& v, auto const& vel, auto const& p) {
        v = coef1;
        v *= p;
        v += vel;
    };
    setvel(v0, velx2, pEx);
    setvel(v1, vely2, pEy);
    setvel(v2, velz2, pEz);
}

template<typename Particles_t, typename... Args>
auto avx_ad_acc_ad(Particles_t& ps, Args&... args) _PHARE_ALL_FN_
{
    auto const& [em, interpolator, layout, halfDtOverDl, dto2m] = std::forward_as_tuple(args...);

    using Float                    = std::decay_t<decltype(halfDtOverDl)>::value_type;
    auto constexpr static simdSize = mkn::avx::Span<Float>::N;

    auto const& siz = ps.size();

    std::size_t i = 0;
    for (; i < siz; i += simdSize)
    {
        // advance_avx(ps, i, halfDtOverDl);
        for (std::size_t j = 0; j < simdSize; ++j)
            advance_noavx(ps, i + j, halfDtOverDl);

        accelerate_avx(ps, i, em, interpolator, layout, halfDtOverDl, dto2m);

        // advance_avx(ps, i, halfDtOverDl);
        for (std::size_t j = 0; j < simdSize; ++j)
            advance_noavx(ps, i + j, halfDtOverDl);
    }

    // do rest
    for (; i < siz; ++i)
    {
        advance_noavx(ps, i, halfDtOverDl);
        accelerate_noavx(ps, i, em, interpolator, layout, halfDtOverDl, dto2m);
        advance_noavx(ps, i, halfDtOverDl);
    }
}

#else

template<typename Particles_t, typename... Args>
auto avx_ad_acc_ad(Particles_t& ps, Args&... args) _PHARE_ALL_FN_
{
    auto const& [em, interpolator, layout, halfDtOverDl, dto2m] = std::forward_as_tuple(args...);

    auto const& siz = ps.size();
    for (std::size_t i = 0; i < siz; ++i)
    {
        advance_noavx(ps, i, halfDtOverDl);
        accelerate_noavx(ps, i, em, interpolator, layout, halfDtOverDl, dto2m);
        advance_noavx(ps, i, halfDtOverDl);
    }
}


#endif // PHARE_HAVE_MKN_AVX


} // namespace PHARE::core::boris


#endif /*PHARE_CORE_PUSHER_BORIS_BASICS_HPP*/
