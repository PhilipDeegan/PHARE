#ifndef PHARE_FLUID_PARTICLE_INITIALIZER_H
#define PHARE_FLUID_PARTICLE_INITIALIZER_H

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/hybrid/hybrid_quantities.h"
#include "core/utilities/types.h"
#include "core/data/ions/particle_initializers/particle_initializer.h"
#include "core/data/particles/particle.h"
#include "initializer/data_provider.h"
#include "core/utilities/point/point.h"


namespace PHARE
{
namespace core
{
    void maxwellianVelocity(std::array<double, 3> V, std::array<double, 3> Vth,
                            std::mt19937_64& generator, std::array<double, 3>& partVelocity);


    std::array<double, 3> basisTransform(std::array<std::array<double, 3>, 3> const& basis,
                                         std::array<double, 3> const& vec);

    void localMagneticBasis(std::array<double, 3> B, std::array<std::array<double, 3>, 3>& basis);


    /** @brief a MaxwellianParticleInitializer is a ParticleInitializer that loads particles from a
     * local Maxwellian distribution given density, bulk velocity and thermal velocity profiles.
     */
    template<typename ParticleArray, typename GridLayout>
    class MaxwellianParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
    {
    public:
        static constexpr auto dimension = GridLayout::dimension;
        using InputFunction             = initializer::InitFunction<dimension>;

        MaxwellianParticleInitializer(InputFunction density,
                                      std::array<InputFunction, 3> bulkVelocity,
                                      std::array<InputFunction, 3> thermalVelocity,
                                      double particleCharge, std::uint32_t nbrParticlesPerCell,
                                      std::optional<std::size_t> seed = {},
                                      Basis basis                     = Basis::Cartesian,
                                      std::array<InputFunction, 3> magneticField
                                      = {nullptr, nullptr, nullptr})
            : density_{density}
            , bulkVelocity_{bulkVelocity}
            , thermalVelocity_{thermalVelocity}
            , magneticField_{magneticField}
            , particleCharge_{particleCharge}
            , nbrParticlePerCell_{nbrParticlesPerCell}
            , basis_{basis}
            , rngSeed_{seed}
        {
        }


        /**
         * @brief load particles in a ParticleArray in a domain defined by the given layout
         */
        void loadParticles(ParticleArray& particles, GridLayout const& layout) const override
        {
            loadParticlesND(density_, particles, layout);
        }


        virtual ~MaxwellianParticleInitializer() = default;


        static std::mt19937_64 getRNG(std::optional<std::size_t> const& seed)
        {
            if (!seed.has_value())
            {
                std::random_device randSeed;
                std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                                       randSeed(), randSeed(), randSeed(), randSeed()};
                return std::mt19937_64(seed_seq);
            }
            return std::mt19937_64(*seed);
        }

    private:
        template<typename Particles /*unnecessary but prevents unused functions from existing*/>
        inline void loadParticlesND(initializer::InitFunction<1> const&, Particles& particles,
                                    GridLayout const& layout) const;
        template<typename Particles>
        inline void loadParticlesND(initializer::InitFunction<2> const&, Particles& particles,
                                    GridLayout const& layout) const;
        template<typename Particles>
        inline void loadParticlesND(initializer::InitFunction<3> const&, Particles& particles,
                                    GridLayout const& layout) const;


        InputFunction density_;
        std::array<InputFunction, 3> bulkVelocity_;
        std::array<InputFunction, 3> thermalVelocity_;
        std::array<InputFunction, 3> magneticField_;

        double particleCharge_;
        std::uint32_t nbrParticlePerCell_;
        Basis basis_;
        std::optional<std::size_t> rngSeed_;
    };
} // namespace core
} // namespace PHARE


namespace PHARE::core
{
namespace
{ // not strictly maxwellian but some what generic
    class MaxwellianInitFunctions
    {
    public:
        template<typename Function, typename FunctionArray, typename Basis, typename... Args>
        MaxwellianInitFunctions(Function& density, FunctionArray& bulkVelocity,
                                FunctionArray& thermalVelocity, FunctionArray& magneticField,
                                Basis const& basis, Args const&... args)
            : _n{density(args...)}
        {
            for (std::uint32_t i = 0; i < 3; i++)
            {
                _V[i]   = bulkVelocity[i](args...);
                _Vth[i] = thermalVelocity[i](args...);
            }
            if (basis == Basis::Magnetic)
                for (std::uint32_t i = 0; i < 3; i++)
                    _B[i] = magneticField[i](args...);
        }

        std::array<double const*, 3> B() const { return ptrs(_B); }

        auto operator()() const { return std::make_tuple(_n->data(), ptrs(_V), ptrs(_Vth)); }

    private:
        std::array<double const*, 3>
        ptrs(std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> const& v) const
        {
            return {v[0]->data(), v[1]->data(), v[2]->data()};
        }
        std::shared_ptr<PHARE::core::Span<double>> const _n;
        std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> _B, _V, _Vth;
    };
} // namespace


template<typename ParticleArray, typename GridLayout>
template<typename Particles>
void MaxwellianParticleInitializer<ParticleArray, GridLayout>::loadParticlesND(
    initializer::InitFunction<1> const&, Particles& particles, GridLayout const& layout) const
{
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);
    auto const [x]        = layout.domainGrids(QtyCentering::primal, [&](auto const&... args) {
        return layout.cellCenteredCoordinates(args...);
    });

    assert(x.size() == (ix1 - ix0));

    particles.reserve(x.size() * nbrParticlePerCell_);

    MaxwellianInitFunctions const fns{density_,       bulkVelocity_, thermalVelocity_,
                                      magneticField_, basis_,        x};

    auto const [n, V, Vth] = fns();
    auto const cellVolume  = layout.cellVolume();
    auto generator         = getRNG(rngSeed_);
    std::size_t cell_idx   = 0;

    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        auto const cellWeight   = n[cell_idx] * cellVolume / nbrParticlePerCell_;
        auto const AMRCellIndex = layout.localToAMR(Point{ix});
        ParticleDeltaDistribution randPosX;
        std::array<double, 3> particleVelocity;
        std::array<std::array<double, 3>, 3> basis;

        if (basis_ == Basis::Magnetic)
        {
            auto const B = fns.B();
            localMagneticBasis({B[0][cell_idx], B[1][cell_idx], B[2][cell_idx]}, basis);
        }

        for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
        {
            maxwellianVelocity({V[0][cell_idx], V[1][cell_idx], V[2][cell_idx]},
                               {Vth[0][cell_idx], Vth[1][cell_idx], Vth[2][cell_idx]}, //
                               generator, particleVelocity);

            if (basis_ == Basis::Magnetic)
                particleVelocity = basisTransform(basis, particleVelocity);

            particles.emplace_back(
                typename Particles::value_type{cellWeight,
                                               particleCharge_,
                                               AMRCellIndex.template toArray<int>(),
                                               {{randPosX(generator)}},
                                               particleVelocity});
        }
        cell_idx++;
    }
}


template<typename ParticleArray, typename GridLayout>
template<typename Particles>
void MaxwellianParticleInitializer<ParticleArray, GridLayout>::loadParticlesND(
    initializer::InitFunction<2> const&, Particles& particles, GridLayout const& layout) const
{
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);
    auto const [iy0, iy1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Y);
    auto const [x, y]     = layout.domainGrids(QtyCentering::primal, [&](auto const&... args) {
        return layout.cellCenteredCoordinates(args...);
    });

    assert(x.size() == y.size());
    particles.reserve(x.size() * nbrParticlePerCell_);

    MaxwellianInitFunctions const fns{
        density_, bulkVelocity_, thermalVelocity_, magneticField_, basis_, x, y};

    auto const [n, V, Vth] = fns();
    auto const cellVolume  = layout.cellVolume();
    auto generator         = getRNG(rngSeed_);
    std::size_t cell_idx   = 0;

    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        for (std::uint32_t iy = iy0; iy < iy1; ++iy)
        {
            auto const cellWeight   = n[cell_idx] * cellVolume / nbrParticlePerCell_;
            auto const AMRCellIndex = layout.localToAMR(Point{ix, iy});
            ParticleDeltaDistribution randPosX, randPosY;
            std::array<double, 3> particleVelocity;
            std::array<std::array<double, 3>, 3> basis;

            if (basis_ == Basis::Magnetic)
            {
                auto const B = fns.B();
                localMagneticBasis({B[0][cell_idx], B[1][cell_idx], B[2][cell_idx]}, basis);
            }

            for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
            {
                maxwellianVelocity({V[0][cell_idx], V[1][cell_idx], V[2][cell_idx]},       //
                                   {Vth[0][cell_idx], Vth[1][cell_idx], Vth[2][cell_idx]}, //
                                   generator, particleVelocity);

                if (basis_ == Basis::Magnetic)
                    particleVelocity = basisTransform(basis, particleVelocity);

                particles.emplace_back(
                    typename Particles::value_type{cellWeight,
                                                   particleCharge_,
                                                   AMRCellIndex.template toArray<int>(),
                                                   {{randPosX(generator), randPosY(generator)}},
                                                   particleVelocity});
            }
            cell_idx++;
        }
    }
}

template<typename ParticleArray, typename GridLayout>
template<typename Particles>
void MaxwellianParticleInitializer<ParticleArray, GridLayout>::loadParticlesND(
    initializer::InitFunction<3> const&, Particles& particles, GridLayout const& layout) const
{
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);
    auto const [iy0, iy1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Y);
    auto const [iz0, iz1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Z);
    auto const [x, y, z]  = layout.domainGrids(QtyCentering::primal, [&](auto const&... args) {
        return layout.cellCenteredCoordinates(args...);
    });

    assert(x.size() == y.size() and y.size() == z.size());
    particles.reserve(x.size() * nbrParticlePerCell_);

    MaxwellianInitFunctions const fns{
        density_, bulkVelocity_, thermalVelocity_, magneticField_, basis_, x, y, z};

    auto const [n, V, Vth] = fns();
    auto const cellVolume  = layout.cellVolume();
    auto generator         = getRNG(rngSeed_);
    std::size_t cell_idx   = 0;

    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        for (std::uint32_t iy = iy0; iy < iy1; ++iy)
        {
            for (std::uint32_t iz = iz0; iz < iz1; ++iz)
            {
                auto const cellWeight   = n[cell_idx] * cellVolume / nbrParticlePerCell_;
                auto const AMRCellIndex = layout.localToAMR(Point{ix, iy, iz});
                ParticleDeltaDistribution randPosX, randPosY, randPosZ;
                std::array<double, 3> particleVelocity;
                std::array<std::array<double, 3>, 3> basis;

                if (basis_ == Basis::Magnetic)
                {
                    auto const B = fns.B();
                    localMagneticBasis({B[0][cell_idx], B[1][cell_idx], B[2][cell_idx]}, basis);
                }

                for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
                {
                    maxwellianVelocity({V[0][cell_idx], V[1][cell_idx], V[2][cell_idx]},       //
                                       {Vth[0][cell_idx], Vth[1][cell_idx], Vth[2][cell_idx]}, //
                                       generator, particleVelocity);

                    if (basis_ == Basis::Magnetic)
                        particleVelocity = basisTransform(basis, particleVelocity);

                    particles.emplace_back(typename Particles::value_type{
                        cellWeight,
                        particleCharge_,
                        AMRCellIndex.template toArray<int>(),
                        {{randPosX(generator), randPosY(generator), randPosZ(generator)}},
                        particleVelocity});
                }
                cell_idx++;
            }
        }
    }
}


} // namespace PHARE::core


#endif
