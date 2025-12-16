#ifndef PHARE_CORE_PUSHER_PORIS_HPP
#define PHARE_CORE_PUSHER_PORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "boris.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>

namespace PHARE::core
{




template<std::size_t dim, typename ParticleArray_t, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class PorisPusher
{
    using Particle_t    = ParticleArray_t::value_type;
    using OnChangeIcell = std::function<bool(ParticleArray_t&, Particle_t&,
                                             std::array<int, dim> const&, std::size_t const)>;


public:
    template<bool _copy_particle = true, typename Population, typename Boxing_t>
    void move_level_ghost(Population& pop, Electromag const& em, Interpolator& interpolator,
                          Boxing_t const& boxing)
    {
        PHARE_LOG_SCOPE(3, "Poris::move_level_ghost");

        using ParticleLoop_t = std::conditional_t<_copy_particle, Particle_t, Particle_t&>;

        auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
        auto const& layout            = boxing.layout;
        double const dto2m            = 0.5 * dt_ / pop.mass();

        auto& domain      = pop.domainParticles();
        auto& level_ghost = pop.levelGhostParticles();

        for (std::size_t i = 0; i < level_ghost.size(); ++i) // size might change on iteration!
        {
            ParticleLoop_t particle = level_ghost[i];
            auto const oldCell      = particle.iCell;

            particle.iCell = boris::advance(particle, halfDtOverDl_);
            boris::accelerate(particle, interpolator(particle, em, layout), dto2m);
            particle.iCell = boris::advance(particle, halfDtOverDl_);

            if (oldCell == particle.iCell)
                continue;

            bool const isInDomainBox      = boxing.isInDomainBox(particle);
            bool const should_interpolote = isInDomainBox;

            if (should_interpolote)
                interpolator.particleToMesh( //
                    particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

            if constexpr (!_copy_particle)
            {
                if (isInDomainBox)
                {
                    domain.push_back(particle);
                    level_ghost.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                    continue;
                }

                bool const isInNonLevelGhostBox
                    = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
                bool const isInLevelGhostBox = !isInNonLevelGhostBox;

                if (isInLevelGhostBox)
                    level_ghost.change_icell(particle, oldCell, i);
                else
                {
                    level_ghost.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                }
            }
        }
    }

    template<bool _copy_particle = true, typename Population, typename Boxing_t>
    void move_domain(Population& pop, Electromag const& em, Interpolator& interpolator,
                     Boxing_t const& boxing)
    {
        PHARE_LOG_SCOPE(3, "Poris::move_domain");

        using ParticleLoop_t = std::conditional_t<_copy_particle, Particle_t, Particle_t&>;

        auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();
        auto const& layout            = boxing.layout;
        double const dto2m            = 0.5 * dt_ / pop.mass();

        auto& domain      = pop.domainParticles();
        auto& patch_ghost = pop.patchGhostParticles();

        for (std::size_t i = 0; i < domain.size(); ++i) // size might change on iteration!
        {
            ParticleLoop_t particle = domain[i];
            auto const oldCell      = particle.iCell;

            particle.iCell = boris::advance(particle, halfDtOverDl_);
            boris::accelerate(particle, interpolator(particle, em, layout), dto2m);
            particle.iCell = boris::advance(particle, halfDtOverDl_);

            if (oldCell == particle.iCell)
            {
                interpolator.particleToMesh( //
                    particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

                continue;
            }

            bool const isInDomainBox = boxing.isInDomainBox(particle);
            bool const isInNonLevelGhostBox
                = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
            bool const should_interpolote = isInNonLevelGhostBox;

            if (!should_interpolote)
            {
                if constexpr (!_copy_particle)
                {
                    domain.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                }
                continue;
            }

            interpolator.particleToMesh( //
                particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

            if constexpr (!_copy_particle)
                domain.change_icell(particle, oldCell, i);
        }

        // move to patch_ghost
        if constexpr (!_copy_particle)
        {
            //
            auto range = makeIndexRange(domain);
            range      = domain.partition(
                range, [&](auto const& cell) { return core::isIn(cell, boxing.domainBox); });

            auto const not_in_domain = makeRange(domain, range.iend(), domain.size());
            patch_ghost.reserve(patch_ghost.size() + not_in_domain.size());
            std::copy(not_in_domain.begin(), not_in_domain.end(), std::back_inserter(patch_ghost));
            domain.erase(not_in_domain);
        }

        for (auto const& p : domain)
            if (!isIn(p, boxing.domainBox))
                throw std::runtime_error("invalid domain");
        for (auto const& p : patch_ghost)
            if (isIn(p, boxing.domainBox))
                throw std::runtime_error("invalid patch ghost");
    }

    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts)
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }


private:
    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
