#ifndef PHARE_INO_UPDATER_HPP
#define PHARE_INO_UPDATER_HPP


#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/pusher/poris.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/moments/moments.hpp"


#include "ion_updater.hpp"



namespace PHARE::core
{


template<typename Ions, typename Electromag, typename GridLayout>
class InoUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box           = PHARE::core::Box<int, dimension>;
    using Interpolator  = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField      = Ions::vecfield_type;
    using ParticleArray = Ions::particle_array_type;
    using Particle_t    = ParticleArray::Particle_t;
    using PartIterator  = ParticleArray::iterator;

    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher            = PorisPusher<dimension, ParticleArray, Electromag, Interpolator,
                                          BoundaryCondition, GridLayout>;

private:
    Pusher pusher_;
    Interpolator interpolator_;

public:
    InoUpdater() = default;
    InoUpdater(auto const& /*dict*/) {}

    template<typename Boxing_t>
    void updatePopulations(Ions& ions, Electromag const& em, Boxing_t const& boxing, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions);


    void reset() { /* noop */ }


private:
    template<typename Boxing_t>
    void updateCopy(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Boxing_t>
    void updateInplace(Ions& ions, Electromag const& em, Boxing_t const& boxing);
};



template<typename Ions, typename Electromag, typename GridLayout>
void InoUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void InoUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 Boxing_t const& boxing, double dt,
                                                                 UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "InoUpdater::updatePopulations");

    resetMoments(ions);
    pusher_.setMeshAndTimeStep(boxing.layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
        updateCopy(ions, em, boxing);

    else
        updateInplace(ions, em, boxing);
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void InoUpdater<Ions, Electromag, GridLayout>::updateCopy(Ions& ions, Electromag const& em,
                                                          Boxing_t const& boxing)
{
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    bool constexpr copy_particle  = true;

    PHARE_LOG_SCOPE(3, "InoUpdater::updateCopy");

    for (auto& pop : ions)
    {
        pusher_.template move_domain<copy_particle>(pop, em, interpolator_, boxing);
        pusher_.template move_level_ghost<copy_particle>(pop, em, interpolator_, boxing);
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void InoUpdater<Ions, Electromag, GridLayout>::updateInplace(Ions& ions, Electromag const& em,
                                                             Boxing_t const& boxing)
{
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
    bool constexpr copy_particle  = false;

    PHARE_LOG_SCOPE(3, "InoUpdater::updateInplace");

    for (auto& pop : ions)
    {
        pusher_.template move_domain<copy_particle>(pop, em, interpolator_, boxing);
        pusher_.template move_level_ghost<copy_particle>(pop, em, interpolator_, boxing);
    }
}



} // namespace PHARE::core

#endif // PHARE_INO_UPDATER_HPP
