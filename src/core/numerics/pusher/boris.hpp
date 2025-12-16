#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/pusher/pusher.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <exception>


namespace PHARE::core::boris
{


struct MoveTwoCellException : std::exception
{
    MoveTwoCellException(double const d, double const v)
        : delta{d}
        , vel{v}
    {
    }

    double delta, vel;
};


template<typename Particle_t, typename Float, std::size_t dim>
auto static advance(Particle_t& p, std::array<Float, dim> const& halfDtOverDl)
{
    std::array<int, dim> newCell;
    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        double const delta = p.delta[iDim] + static_cast<double>(halfDtOverDl[iDim] * p.v[iDim]);

        if (std::abs(delta) > 2)
            throw MoveTwoCellException{delta, p.v[iDim]};

        auto const iCell = static_cast<int>(std::floor(delta));

        p.delta[iDim] = delta - iCell;
        newCell[iDim] = iCell + p.iCell[iDim];
    }
    return newCell;
}



/** Accelerate the particles in rangeIn and put the new velocity in rangeOut
 */
template<typename Particle, typename ParticleEB, typename Float>
void accelerate(Particle& p, ParticleEB const& eb, Float const& dto2m)
{
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    Float const coef1 = p.charge * dto2m;

    // We now apply the 3 steps of the BORIS PUSHER
    // 1st half push of the electric field
    Float const velx1 = p.v[0] + coef1 * pEx;
    Float const vely1 = p.v[1] + coef1 * pEy;
    Float const velz1 = p.v[2] + coef1 * pEz;

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
    p.v[0] = velx2 + coef1 * pEx;
    p.v[1] = vely2 + coef1 * pEy;
    p.v[2] = velz2 + coef1 * pEz;
}


} // namespace PHARE::core::boris


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
        PHARE_LOG_SCOPE(3, "Boris::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        prePushStep_(rangeIn, rangeOut);


        rangeOut = firstSelector(rangeOut);

        double const dto2m = 0.5 * dt_ / mass;
        for (auto idx = rangeOut.ibegin(); idx < rangeOut.iend(); ++idx)
        {
            auto& currPart = rangeOut.array()[idx];

            //  get electromagnetic fields interpolated on the particles of rangeOut stop at newEnd.
            //  get the particle velocity from t=n to t=n+1
            auto const& local_em = interpolator(currPart, emFields, layout);
            boris::accelerate(currPart, local_em, dto2m);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            try
            {
                postPushStep_(rangeOut, idx);
            }
            catch (DictionaryException const& bex)
            {
                auto ex            = bex;
                auto const& [e, b] = local_em;
                for (std::uint16_t i = 0; i < 3; ++i)
                    ex("E_" + std::to_string(i), std::to_string(e[i]));
                for (std::uint16_t i = 0; i < 3; ++i)
                    ex("B_" + std::to_string(i), std::to_string(b[i]));
                ex("level", std::to_string(layout.levelNumber()));
                throw ex;
            }
        }

        return secondSelector(rangeOut);
    }



    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename Particle>
    auto advancePosition_(Particle const& partIn, Particle& partOut)
    {
        std::array<int, dim> newCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double const delta
                = partIn.delta[iDim] + static_cast<double>(halfDtOverDl_[iDim] * partIn.v[iDim]);

            if (std::abs(delta) > 2)
                throw boris::MoveTwoCellException{delta, partIn.v[iDim]};

            auto const iCell = static_cast<int>(std::floor(delta));

            partOut.delta[iDim] = delta - iCell;
            newCell[iDim]       = iCell + partIn.iCell[iDim];
        }
        return newCell;
    }


    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the ParticleSelector
     */
    void prePushStep_(ParticleRange const& rangeIn, ParticleRange& rangeOut)
    {
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

            outParticles[outIdx].charge = inParticles[inIdx].charge;
            outParticles[outIdx].weight = inParticles[inIdx].weight;
            outParticles[outIdx].v      = inParticles[inIdx].v;

            try
            {
                auto newCell = advancePosition_(inParticles[inIdx], outParticles[outIdx]);
                if (newCell != inParticles[inIdx].iCell)
                    outParticles.change_icell(newCell, outIdx);
            }
            catch (boris::MoveTwoCellException const& e)
            {
                std::stringstream ss;
                ss << "PrePush Particle moved 2 cells with delta/vel: ";
                ss << e.delta << "/" << e.vel << std::endl;
                DictionaryException ex{"cause", ss.str()};
                throw ex;
            }
        }
    }

    void postPushStep_(ParticleRange& range, std::size_t idx)
    {
        try
        {
            auto& particles = range.array();
            auto newCell    = advancePosition_(particles[idx], particles[idx]);
            if (newCell != particles[idx].iCell)
                particles.change_icell(newCell, idx);
        }
        catch (boris::MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PostPush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            throw DictionaryException{}("cause", ss.str());
        }
    }



    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
