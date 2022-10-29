#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include "core/numerics/pusher/pusher.hpp"
#include "core/utilities/range/range.hpp"
#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/data/particles/particle.hpp"

namespace PHARE::core
{
template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
    : public Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>;

private:
    using ParticleSelector = typename Super::ParticleSelector;

public:
    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    /** see Pusher::move() documentation*/
#if 0
    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double mass, Interpolator& interpolator,
                       ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                       GridLayout const& layout) override
    {
            // push the particles of half a step
            // rangeIn : t=n, rangeOut : t=n+1/Z
            // get a pointer on the first particle of rangeOut that leaves the patch
            auto firstLeaving
                = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

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
#endif


    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double mass, Interpolator& interpolator,
                       GridLayout const& layout, ParticleSelector firstSelector,
                       ParticleSelector secondSelector) override
    {
        PHARE_LOG_SCOPE("Boris::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        pushStep_(rangeIn, rangeOut, PushStep::PrePush);

        rangeOut = firstSelector(rangeOut);

        //  get electromagnetic fields interpolated on the particles of rangeOut
        //  stop at newEnd.
        interpolator.meshToParticle(rangeOut, emFields, layout);

        //  get the particle velocity from t=n to t=n+1
        accelerate_(rangeOut, rangeOut, mass);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        pushStep_(rangeOut, rangeOut, PushStep::PostPush);

        return secondSelector(rangeOut);
    }




    /** see Pusher::move() documentation*/
    virtual void setMeshAndTimeStep(std::array<double, dim> ms, double ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    enum class PushStep { PrePush, PostPush };

    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename Particle>
    auto advancePosition_(Particle const& partIn, Particle& partOut)
    {
        std::array<int, dim> newCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double delta = partIn.delta()[iDim]
                           + static_cast<double>(halfDtOverDl_[iDim] * partIn.v()[iDim]);

            double iCell = std::floor(delta);
            if (std::abs(delta) > 2)
            {
                PHARE_LOG_ERROR("Error, particle moves more than 1 cell, delta >2");
            }
            partOut.delta()[iDim] = delta - iCell;
            newCell[iDim]         = static_cast<int>(iCell + partIn.iCell()[iDim]);
        }
        return newCell;
    }



    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the ParticleSelector
     */
    void pushStep_(ParticleRange const& rangeIn, ParticleRange& rangeOut, PushStep step)
    {
        using ParticleArray = std::decay_t<decltype(rangeIn.array())>;

        auto& inParticles  = rangeIn.array();
        auto& outParticles = rangeOut.array();
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
             ++inIdx, ++outIdx)
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
                outParticles[outIdx].charge() = inParticles[inIdx].charge();
                outParticles[outIdx].weight() = inParticles[inIdx].weight();
                outParticles[outIdx].v()      = inParticles[inIdx].v();
            }
            auto newCell = advancePosition_(inParticles[inIdx], outParticles[outIdx]);

            if constexpr (ParticleArray::is_mapped)
                if (newCell != inParticles[inIdx].iCell())
                    outParticles.change_icell(newCell, outIdx);
        }
    }



    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename ParticleRangeIn, typename ParticleRangeOut>
    void accelerate_(ParticleRangeIn const& inputParticles, ParticleRangeOut& outputParticles,
                     double const mass)
    {
        double dto2m = 0.5 * dt_ / mass;

        auto out_idx        = outputParticles.begin().idx();
        auto& out_particles = outputParticles.begin()();

        auto in_start      = inputParticles.begin().idx();
        auto in_end        = inputParticles.end().idx();
        auto& in_particles = inputParticles.begin()();

        for (auto in_idx = in_start; in_idx < in_end; ++in_idx)
        {
            auto const& E = in_particles.E(in_idx);
            auto const& B = in_particles.B(in_idx);

            auto const& [Ex, Ey, Ez] = E;
            auto const& [Bx, By, Bz] = B;

            double coef1 = in_particles.charge(in_idx) * dto2m;

            // We now apply the 3 steps of the BORIS PUSHER

            // 1st half push of the electric field
            auto const& in_v                = in_particles.v(in_idx);
            std::array<double, 3> const vel = {in_v[0] + coef1 * Ex, //
                                               in_v[1] + coef1 * Ey, //
                                               in_v[2] + coef1 * Ez};

            // preparing variables for magnetic rotation
            double const rx = coef1 * Bx;
            double const ry = coef1 * By;
            double const rz = coef1 * Bz;

            double const rx2  = rx * rx;
            double const ry2  = ry * ry;
            double const rz2  = rz * rz;
            double const rxry = rx * ry;
            double const rxrz = rx * rz;
            double const ryrz = ry * rz;

            double const invDet = 1. / (1. + rx2 + ry2 + rz2);

            // preparing rotation matrix due to the magnetic field
            // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
            double const mxx = 1. + rx2 - ry2 - rz2;
            double const mxy = 2. * (rxry + rz);
            double const mxz = 2. * (rxrz - ry);

            double const myx = 2. * (rxry - rz);
            double const myy = 1. + ry2 - rx2 - rz2;
            double const myz = 2. * (ryrz + rx);

            double const mzx = 2. * (rxrz + ry);
            double const mzy = 2. * (ryrz - rx);
            double const mzz = 1. + rz2 - rx2 - ry2;

            // magnetic rotation
            double const velx2 = (mxx * vel[0] + mxy * vel[1] + mxz * vel[2]) * invDet;
            double const vely2 = (myx * vel[0] + myy * vel[1] + myz * vel[2]) * invDet;
            double const velz2 = (mzx * vel[0] + mzy * vel[1] + mzz * vel[2]) * invDet;

            // 2nd half push of the electric field
            // Update particle velocity
            out_particles.v(out_idx) = {velx2 + coef1 * Ex, //
                                        vely2 + coef1 * Ey, //
                                        velz2 + coef1 * Ez};
            ++out_idx;
        }
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
