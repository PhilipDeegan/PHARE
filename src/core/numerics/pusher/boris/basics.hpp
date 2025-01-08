#ifndef PHARE_CORE_PUSHER_BORIS_BASICS_HPP
#define PHARE_CORE_PUSHER_BORIS_BASICS_HPP


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



} // namespace PHARE::core::boris


#endif /*PHARE_CORE_PUSHER_BORIS_BASICS_HPP*/
