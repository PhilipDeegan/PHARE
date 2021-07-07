#ifndef PHARE_CORE_PUSHER_BORIS_H
#define PHARE_CORE_PUSHER_BORIS_H

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "core/data/particles/contiguous.h"
#include "core/numerics/pusher/pusher.h"
#include "core/utilities/range/range.h"
#include "core/errors.h"
#include "core/logger.h"


template<typename Float>
void c_boris_advancePosition(size_t dim, Float const* vIn, Float const* deltaIn, int const* iCellIn,
                             Float* deltaOut, int* iCellOut, Float const* halfDtOverDl)
{
    // push the particle
    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        Float delta = deltaIn[iDim] + static_cast<Float>(halfDtOverDl[iDim] * vIn[iDim]);
        Float iCell = std::floor(delta);
        if (std::abs(delta) > 2)
        {
            throw std::runtime_error("Error, particle moves more than 1 cell, delta >2");
        }
        deltaOut[iDim] = delta - iCell;
        iCellOut[iDim] = static_cast<int>(iCell + iCellIn[iDim]);
    }
}

template<typename Float>
void c_boris_accelerate(Float const& chargeIn, Float const* vIn, Float const* EIn,
                        Float const* const BIn, Float* vOut, Float dto2m)
{
    constexpr Float one = 1;
    constexpr Float two = 2;

    Float coef1 = chargeIn * dto2m;

    // We now apply the 3 steps of the BORIS PUSHER

    // 1st half push of the electric field
    Float velx1 = vIn[0] + coef1 * EIn[0]; // see Fused Multiple/Add (FMA)
    Float vely1 = vIn[1] + coef1 * EIn[1];
    Float velz1 = vIn[2] + coef1 * EIn[2];


    // preparing variables for magnetic rotation
    Float const rx = coef1 * BIn[0];
    Float const ry = coef1 * BIn[1];
    Float const rz = coef1 * BIn[2];

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


    // 2nd half push of the electric field
    // Update particle velocity
    vOut[0] = velx2 + coef1 * EIn[0];
    vOut[1] = vely2 + coef1 * EIn[1];
    vOut[2] = velz2 + coef1 * EIn[2];
}


namespace PHARE::core
{
template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
    : public Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>;
    using Float            = typename GridLayout::Float;
    using ParticleSelector = typename Super::ParticleSelector;
    using ParticleRange    = Range<ParticleIterator>;

    static constexpr Float half = .5;


    /** see Pusher::move() domentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, Float mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                          GridLayout const& layout) override
    {
        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/Z
        // get a pointer on the first particle of rangeOut that leaves the patch
        auto firstLeaving = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

        // apply boundary condition on the particles in [firstLeaving, rangeOut.end[
        // that actually leave through a physical boundary condition
        // get a pointer on the new end of rangeOut. Particles passed newEnd
        // are those that have left the patch through a non-physical boundary
        // they should be discarded now
        auto newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

        rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

        // get electromagnetic fields interpolated on the particles of rangeOut
        // stop at newEnd.
        interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

        // get the particle velocity from t=n to t=n+1
        accelerate_(rangeOut, rangeOut, mass);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

        // apply BC on the leaving particles that leave through physical BC
        // and get pointer on new End, discarding particles leaving elsewhere
        newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

        rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

        return rangeOut.end();
    }


    /** see Pusher::move() domentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, Float mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving,
                          GridLayout const& layout) override
    {
        PHARE_LOG_SCOPE("Boris::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        pushStep_(rangeIn, rangeOut, PushStep::PrePush);

        // get electromagnetic fields interpolated on the particles of rangeOut
        // stop at newEnd.
        interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

        // get the particle velocity from t=n to t=n+1
        accelerate_(rangeOut, rangeOut, mass);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        auto firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

        rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

        return rangeOut.end();
    }


    /** see Pusher::move() domentation*/
    virtual void setMeshAndTimeStep(std::array<Float, dim> ms, Float ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](Float& x) { return half * ts / x; });
        dt_ = ts;
    }



private:
    enum class PushStep { PrePush, PostPush };

    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename ParticleIter>
    void advancePosition_(ParticleIter const& partIn, ParticleIter& partOut)
    {
        // push the particle
        c_boris_advancePosition(dim, partIn.v.data(), partIn.delta.data(), partIn.iCell.data(),
                                partOut.delta.data(), partOut.iCell.data(), halfDtOverDl_.data());
    }



    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the ParticleSelector
     */
    template<typename ParticleRangeIn, typename ParticleRangeOut>
    void pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut, PushStep step)
    {
        auto currentOut = rangeOut.begin();
        for (auto& currentIn : rangeIn)
        {
            // in the first push, this is the first time
            // we push to rangeOut, which contains crap
            // the push will only touch the particle position
            // but the next step being the acceleration of
            // rangeOut, we need to copy rangeIn weights, charge
            // and velocity. This is done here although
            // not strictly speaking this function's business
            // to take advantage that we're already looping
            // over rangeIn particles.
            if (step == PushStep::PrePush)
            {
                currentOut->charge = currentIn.charge;
                currentOut->weight = currentIn.weight;
                currentOut->v      = currentIn.v;
            }
            // push the particle
            advancePosition_(currentIn, *currentOut);
            currentOut++;
        }
    }

    template<typename ParticleRangeIn, typename ParticleRangeOut>
    auto pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut,
                   ParticleSelector const& particleIsNotLeaving, PushStep step)
    {
        pushStep_(rangeIn, rangeOut, step);

        // now all particles have been pushed
        // those not satisfying the predicate after the push
        // are found in [newEnd:end[
        // those for which pred is true are in [firstOut,newEnd[
        return std::partition(std::begin(rangeOut), std::end(rangeOut), particleIsNotLeaving);
    }




    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename ParticleRangeIn, typename ParticleRangeOut>
    void accelerate_(ParticleRangeIn inputParticles, ParticleRangeOut outputParticles, Float mass)
    {
        Float dto2m = half * dt_ / mass;

        auto currentOut = outputParticles.begin();

        for (auto const& currentIn : inputParticles)
        {
            std::vector<Float const*> const E = {&currentIn.Ex, &currentIn.Ey, &currentIn.Ez};
            std::vector<Float const*> const B = {&currentIn.Bx, &currentIn.By, &currentIn.Bz};
            c_boris_accelerate(currentIn.charge, currentIn.v.data(), E[0], B[0],
                               currentOut->v.data(), dto2m);

            ++currentOut;
        }
    }


    std::array<Float, dim> halfDtOverDl_;
    Float dt_;
};


} // namespace PHARE::core


#endif
