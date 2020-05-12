#ifndef PHARE_CORE_PUSHER_BORIS_H
#define PHARE_CORE_PUSHER_BORIS_H

#include <array>
#include <cmath>
#include <cstddef>

#include "core/data/particles/contiguous.h"
#include "core/numerics/pusher/pusher.h"
#include "core/utilities/range/range.h"

namespace PHARE::core
{
template<typename Float>
void c_boris_advancePosition(size_t dim, Float const* vIn, float const* deltaIn, int const* iCellIn,
                             float* deltaOut, int* iCellOut, Float const* halfDtOverDl)
{
    // push the particle
    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        float delta = deltaIn[iDim] + static_cast<float>(halfDtOverDl[iDim] * vIn[iDim]);

        float iCell    = std::floor(delta);
        deltaOut[iDim] = delta - iCell;
        iCellOut[iDim] = static_cast<int>(iCell + iCellIn[iDim]);
    }
}

template<typename Float>
// for SIMD this should be split up to a single function per operation - NO LAMBDAS!
void c_boris_accelerate(Float const& chargeIn, Float const* vIn, Float const* EIn,
                        Float const* const BIn, Float* vOut, Float dto2m)
{
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

    Float const invDet = 1. / (1. + rx2 + ry2 + rz2);

    // preparing rotation matrix due to the magnetic field
    // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
    Float const mxx = 1. + rx2 - ry2 - rz2;
    Float const mxy = 2. * (rxry + rz);
    Float const mxz = 2. * (rxrz - ry);

    Float const myx = 2. * (rxry - rz);
    Float const myy = 1. + ry2 - rx2 - rz2;
    Float const myz = 2. * (ryrz + rx);

    Float const mzx = 2. * (rxrz + ry);
    Float const mzy = 2. * (ryrz - rx);
    Float const mzz = 1. + rz2 - rx2 - ry2;

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



template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
         typename ParticleSelector, typename BoundaryCondition, typename GridLayout>
class BorisPusher : public Pusher<dim, ParticleIterator, Electromag, Interpolator, ParticleSelector,
                                  BoundaryCondition, GridLayout>
{
public:
    using Float         = typename GridLayout::float_type;
    using ParticleRange = Range<ParticleIterator>;
    using Particles     = EMContiguousParticles<Float, dim>;

    /** see Pusher::move() domentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, double mass, Interpolator& interpolator,
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
                          Electromag const& emFields, double mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving,
                          GridLayout const& layout) override
    {
        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/Z
        // get a pointer on the first particle of rangeOut that leaves the patch
        auto firstLeaving = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

        rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

        // get electromagnetic fields interpolated on the particles of rangeOut
        // stop at newEnd.
        interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

        // get the particle velocity from t=n to t=n+1
        accelerate_(rangeOut, rangeOut, mass);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

        rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

        return rangeOut.end();
    }


    /** see Pusher::move() domentation*/
    void setMeshAndTimeStep(std::array<double, dim> ms, double ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }

    template<size_t PushStep>
    void move(Particles& particleIn, Particles& particleOut,
              std::vector<std::vector<float>>& elecromag /*Exyz/Bxyz*/)
    {
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
    auto pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut,
                   ParticleSelector const& particleIsNotLeaving, PushStep step)
    {
        auto swapee = rangeOut.end() - 1;
        auto newEnd = rangeOut.end();

        auto currentOut = rangeOut.begin();

        for (auto& currentIn : rangeIn)
        {
            if (step == PushStep::PrePush)
            {
                currentOut->charge = currentIn.charge;
                currentOut->weight = currentIn.weight;
                for (std::size_t i = 0; i < 3; ++i)
                    currentOut->v[i] = currentIn.v[i];
            }
            // push the particle
            advancePosition_(currentIn, *currentOut);

            if (particleIsNotLeaving(*currentOut))
            {
                // we advance the output iterator
                // only if currentOut has not been
                // swapped
                ++currentOut;
            }
            else
            {
                // if the particle satisfies the predicate
                // swap it with the swapee
                // and decrement the swapee

                std::swap(*currentOut, *swapee);
                --newEnd;
                --swapee;
            }
        }

        // now all particles have been pushed
        // those not satisfying the predicate after the push
        // are found in [newEnd:end[
        // those for which pred is true are in [firstOut,newEnd[
        return newEnd;
    }




    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename ParticleRangeIn, typename ParticleRangeOut>
    void accelerate_(ParticleRangeIn inputParticles, ParticleRangeOut outputParticles, double mass)
    {
        double dto2m = 0.5 * dt_ / mass;

        auto currentOut = outputParticles.begin();

        for (auto const& currentIn : inputParticles)
        {
            std::vector<double const*> const E = {&currentIn.Ex, &currentIn.Ey, &currentIn.Ez};
            std::vector<double const*> const B = {&currentIn.Bx, &currentIn.By, &currentIn.Bz};
            c_boris_accelerate(currentIn.charge, currentIn.v.data(), E[0], B[0],
                               currentOut->v.data(), dto2m);

            ++currentOut;
        }
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};


} // namespace PHARE::core


#endif
