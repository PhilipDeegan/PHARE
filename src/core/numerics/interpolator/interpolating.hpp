#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "interpolator.hpp"

namespace PHARE::core
{

template<typename ParticleArray_t, std::size_t interpOrder, bool atomic_ops,
         typename Interpolator_t
         = Interpolator<ParticleArray_t::dimension, interpOrder, atomic_ops>>
class Interpolating
{
    auto constexpr static dim = ParticleArray_t::dimension;

public:
    template<typename Particles, typename VecField, typename GridLayout, typename Field>
    inline void operator()(Particles const& particles, Field& rhoP, Field& rhoC, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
        requires(ParticleArray_t::layout_mode == LayoutMode::AoSMapped)
    {
        Interpolator_t{}(particles, rhoP, rhoC, flux, layout, coef);
    }

    template<typename Particles, typename VecField, typename GridLayout, typename Field>
    inline void operator()(Particles const& particles, Field& rhoP, Field& rhoC, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
        requires(ParticleArray_t::layout_mode != LayoutMode::AoSMapped)
    {
        particleToMesh(particles, layout, rhoP, rhoC, flux, coef);
    }

    // AoSTS: one particle bucket per tile
    template<typename GridLayout, typename VecField, typename Field>
    void particleToMesh(ParticleArray_t const& particles, GridLayout const& layout, Field& rhoP,
                        Field& rhoC, VecField& flux, double coef = 1.)
        requires(ParticleArray_t::layout_mode == LayoutMode::AoSTS)
    {
        using Field_vt    = Field::value_type;
        using VecField_vt = basic::TensorField<Field_vt, 1>;

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto& rhop        = rhoP()[tidx];
            auto& rhoc        = rhoC()[tidx];
            auto F            = flux.template as<VecField_vt>([&](auto& c) { return c()[tidx]; });
            auto const& parts = particles()[tidx];
            for (auto const& p : parts())
                interp_.particleToMesh(p, rhop(), rhoc(), F, rhop.layout(), coef);
        }
    }

    // AoSPCTS: tiles of per-cell buckets
    template<typename GridLayout, typename VecField, typename Field>
    void particleToMesh(ParticleArray_t const& particles, GridLayout const& layout, Field& rhoP,
                        Field& rhoC, VecField& flux, double coef = 1.)
        requires(ParticleArray_t::layout_mode == LayoutMode::AoSPCTS)
    {
        using Field_vt    = Field::value_type;
        using VecField_vt = basic::TensorField<Field_vt, 1>;

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto& rhop   = rhoP()[tidx];
            auto& rhoc   = rhoC()[tidx];
            auto F       = flux.template as<VecField_vt>([&](auto& c) { return c()[tidx]; });
            auto& pctile = particles()[tidx];
            auto& cps    = pctile();
            // full local box (domain + halo), not domain-only: level/patch ghost
            // particles live in the halo and must still be deposited
            for (auto const& bix : cps.local_box())
                for (auto const& p : cps(bix))
                    interp_.particleToMesh(p, rhop(), rhoc(), F, rhop.layout(), coef);
        }
    }

    // AoSPC: flat per-cell buckets, not tiled
    template<typename GridLayout, typename VecField, typename Field>
    void particleToMesh(ParticleArray_t const& particles, GridLayout const& layout, Field& rhoP,
                        Field& rhoC, VecField& flux, double coef = 1.)
        requires(ParticleArray_t::layout_mode == LayoutMode::AoSPC)
    {
        PHARE_LOG_SCOPE(3, "Interpolating::particleToMesh");

        for (auto const& bix : particles.local_box())
            for (auto const& p : particles(bix))
                interp_.particleToMesh(p, rhoP, rhoC, flux, layout, coef);
    }

    // AoS: flat, non-tiled, non-per-cell
    template<typename GridLayout, typename VecField, typename Field>
    void particleToMesh(ParticleArray_t const& particles, GridLayout const& layout, Field& rhoP,
                        Field& rhoC, VecField& flux, double coef = 1.)
        requires(ParticleArray_t::layout_mode == LayoutMode::AoS)
    {
        PHARE_LOG_SCOPE(3, "Interpolating::particleToMesh");

        interp_(particles, rhoP, rhoC, flux, layout, coef);
    }


    Interpolator_t interp_;
};

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
struct MomentumTensorInterpolating
{
    using Interpolator_t = MomentumTensorInterpolator<dim, interpOrder, atomic_ops>;

    Interpolator_t interp_;

public:
    template<typename Particles_t, typename TensorField, typename GridLayout>
    inline void operator()(Particles_t& particles, TensorField& momentumTensor,
                           GridLayout const& layout, double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSMapped)
    {
        assert(momentumTensor.isUsable());
        interp_(particles, momentumTensor, layout, mass);
    }

    // AoSTS: one flat particle bucket per tile
    template<typename Particles_t, typename TensorField, typename GridLayout>
    inline void operator()(Particles_t& particles, TensorField& momentumTensor,
                           GridLayout const& layout, double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSTS)
    {
        assert(momentumTensor.isUsable());

        using Field_vt       = TensorField::field_type::value_type;
        using TensorField_vt = basic::TensorField<Field_vt, 2>;

        assert(particles().size() == momentumTensor[0]().size());

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt
                = momentumTensor.template as<TensorField_vt>([&](auto& c) { return c()[tidx]; });
            assert(mt.isUsable());
            auto const& parts = particles()[tidx];
            interp_(parts(), mt, mt[0].layout(), mass);
        }
    }

    // AoSPCTS: tiles of per-cell buckets - unlike AoSTS a tile has no single flat
    // range, so each cell's bucket is deposited as its own range
    template<typename Particles_t, typename TensorField, typename GridLayout>
    inline void operator()(Particles_t& particles, TensorField& momentumTensor,
                           GridLayout const& layout, double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    {
        assert(momentumTensor.isUsable());

        using Field_vt       = TensorField::field_type::value_type;
        using TensorField_vt = basic::TensorField<Field_vt, 2>;

        assert(particles().size() == momentumTensor[0]().size());

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt
                = momentumTensor.template as<TensorField_vt>([&](auto& c) { return c()[tidx]; });
            assert(mt.isUsable());

            auto& pctile = particles()[tidx];
            auto& cps    = pctile();
            for (auto const& bix : cps.local_box())
                interp_(cps(bix), mt, mt[0].layout(), mass);
        }
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
