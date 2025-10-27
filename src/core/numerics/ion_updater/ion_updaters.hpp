#ifndef PHARE_ION_UPDATERS_HPP
#define PHARE_ION_UPDATERS_HPP


#include "core/numerics/ion_updater/ion_updater.hpp"          // IWYU pragma: keep
#include "core/numerics/ion_updater/ion_updater_multi_ts.hpp" // IWYU pragma: keep
// #include "core/numerics/ion_updater/ion_updater_multi_pc.hpp" // IWYU pragma: keep


namespace PHARE::core
{


template<typename Ions, typename Electromag, typename GridLayout>
struct IonUpdaterImplResolverFns
{
    using ParticleArray_t = Ions::particle_array_type;

    auto static constexpr updater_impl()
    {
        using enum core::LayoutMode;
        if constexpr (any_in(ParticleArray_t::layout_mode, AoSMapped))
            return _as_nullptr_<core::IonUpdater<Ions, Electromag, GridLayout>*>();
#if PHARE_HAVE_MKN_GPU
        else if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS))
            return _as_nullptr_<core::mkn::IonUpdaterMultiTS<Ions, Electromag, GridLayout>*>();
#endif // PHARE_HAVE_MKN_GPU
        else
            throw std::runtime_error("no");
    }

    template<typename T>
    auto static constexpr _as_nullptr_()
    {
        return T{nullptr};
    }
};

template<typename Ions, typename Electromag, typename GridLayout>
struct IonUpdaterImplResolver : IonUpdaterImplResolverFns<Ions, Electromag, GridLayout>
{
    using Super        = IonUpdaterImplResolverFns<Ions, Electromag, GridLayout>;
    using IonUpdater_t = std::decay_t<decltype(*Super::updater_impl())>;
};


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
