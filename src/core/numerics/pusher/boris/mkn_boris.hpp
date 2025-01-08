#ifndef PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP


#if PHARE_HAVE_MKN_GPU

#include "core/utilities/kernels.hpp"


namespace PHARE::core::detail
{
auto static const multi_boris_threads = get_env_as("PHARE_ASYNC_THREADS", std::size_t{5});

} // namespace PHARE::core::detail

namespace PHARE::core
{

// !!naive!!
template<typename ModelViews>
struct MultiBoris
{
    using ModelView    = ModelViews::value_type;
    using GridLayout_t = ModelView::GridLayout_t;

    using ParticleArray_t     = ModelView::ParticleArray_t;
    static constexpr auto dim = ParticleArray_t::dimension;
    using Electromag_t        = ModelView::Electromag_t;
    using Vecfield_t          = Electromag_t::vecfield_type;
    using Field_t             = Vecfield_t::field_type;
    using ParticleArray_v     = typename ParticleArray_t::view_t;
    using Box_t               = Box<std::uint32_t, dim>;
    using Boxes_t             = std::vector<Box_t>;
    using Particles_ptrs      = std::vector<ParticleArray_t*>;
    using StreamLauncher      = gpu::ThreadedBoxStreamLauncher<Boxes_t, Particles_ptrs>;


    static auto _particles(ModelViews& views)
    {
        std::vector<ParticleArray_t*> ptrs;
        auto add = [&](auto& ps) {
            if (ps.size())
                ptrs.emplace_back(&ps);
        };
        auto all = [&](auto&... ps) { (add(ps), ...); };

        for (auto& view : views)
            for (auto& pop : *view.ions)
                all(pop.domainParticles(), pop.patchGhostParticles(), pop.levelGhostParticles());

        return ptrs;
    }

    MultiBoris(double const dt_, ModelViews& _views)
        : dt{dt_}
        , views{_views}
        , particles{_particles(views)}
    {
        for (auto& view : views)
        {
            auto each = [&](auto const& type, auto& ps, auto& pop) mutable {
                if (!ps.size())
                    return;
                particle_type.emplace_back(type);
                auto& domainview = pviews.emplace_back(*ps);
                boxes.emplace_back(domainview.local_box());
                layouts.emplace_back(view.layout);
                ems.emplace_back(*view.em);
                dto2ms.emplace_back(0.5 * dt / pop.mass());
                halfdt.emplace_back(mesh(view.layout.meshSize(), dt));
                rhos.emplace_back(pop.density());
                fluxes.emplace_back(pop.flux());
            };

            for (auto& pop : *view.ions)
            {
                each(0, pop.domainParticles(), pop);
                each(1, pop.patchGhostParticles(), pop);
                each(2, pop.levelGhostParticles(), pop);
            }
        }
    }

    void reset() {}

    double dt;
    ModelViews& views;

    // on host
    Particles_ptrs particles;
    std::vector<std::uint16_t> particle_type;
    Boxes_t boxes;

    // on device
    gpu::Vec_t<ParticleArray_v> pviews;
    gpu::Vec_t<GridLayout_t> layouts;
    gpu::Vec_t<Electromag_t> ems;
    gpu::Vec_t<Field_t> rhos;
    gpu::Vec_t<Vecfield_t> fluxes;
    gpu::Vec_t<std::array<double, dim>> halfdt;
    gpu::Vec_t<double> dto2ms;

    StreamLauncher streamer{particles, boxes, detail::multi_boris_threads};

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};




} // namespace PHARE::core

#endif // PHARE_HAVE_MKN_GPU

#endif /*PHARE_CORE_PUSHER_BORIS_MKN_BORIS_HPP*/
