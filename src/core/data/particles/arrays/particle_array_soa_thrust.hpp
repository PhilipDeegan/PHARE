#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP


#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"


#include <type_traits>


namespace PHARE::core
{


// used in SOA_ParticelArray::operator[](std::size_t) for std algorithm interactions
//  allows copying between vector and array impls without iterators
//   could be a better way maybe but it's not so simple to change
template<std::size_t dim>
using SoAParticle_crt = std::tuple<double const&,                  //  weight
                                   double const&,                  // charge
                                   std::array<int, dim> const&,    // iCell
                                   std::array<double, dim> const&, // delta
                                   std::array<double, 3> const&    // v

                                   >;
template<std::size_t dim>
using SoAParticle_rt = std::tuple<double&,                  //  weight
                                  double&,                  // charge
                                  std::array<int, dim>&,    // iCell
                                  std::array<double, dim>&, // delta
                                  std::array<double, 3>&    // v
                                  >;



} // namespace PHARE::core


namespace PHARE::core::detail
{
template<typename SoAParticles_t>
struct SoAZipParticle;
}


#if PHARE_HAVE_THRUST

#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/iterator/zip_iterator.h>


namespace PHARE::core::detail
{

struct SoAIteratorAdaptor
{
    template<typename Particles>
    auto static make(Particles& particles, std::size_t const& idx = 0) _PHARE_ALL_FN_
    {
        if constexpr (Particles::storage_mode == StorageMode::SPAN)
        {
            return thrust::make_zip_iterator( //
                thrust::make_tuple(           //
                    particles.charge_ + idx,  // 0
                    particles.weight_ + idx,  // 1
                    particles.iCell_ + idx,   // 2
                    particles.delta_ + idx,   // 3
                    particles.v_ + idx        // 4
                    ));
        }

        else
        {
            return thrust::make_zip_iterator(       //
                thrust::make_tuple(                 //
                    particles.charge_.data() + idx, // 0
                    particles.weight_.data() + idx, // 1
                    particles.iCell_.data() + idx,  // 2
                    particles.delta_.data() + idx,  // 3
                    particles.v_.data() + idx       // 4
                    ));
        }
    }
    template<typename T>
    static auto& charge(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<0>(it);
    }
    template<typename T>
    static auto& weight(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<1>(it);
    }

    template<typename T>
    static auto& iCell(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<2>(it);
    }
    template<typename T>
    static auto& delta(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<3>(it);
    }
    template<typename T>
    static auto& v(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<4>(it);
    }
};


template<typename SoAParticles_t>
struct SoAZipParticle
{
    auto constexpr static dim = SoAParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)}
    {
    }

    auto& charge() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& weight() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& iCell() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& delta() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& v() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }

    Iterator it;
};



template<typename SoAParticles_t>
struct SoAZipConstParticle
{
    auto constexpr static dim = SoAParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipConstParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)}
    {
    }

    auto& charge() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }

    Iterator it;
};


struct SoAVXIteratorAdaptor
{
    template<std::size_t dim>
    static auto soaxv_tuple_as_thrust_tuple(auto&& v, std::size_t const idx)
    {
        auto tuple = v.as_tuple();

        auto seq = [&]<std::size_t offset, std::size_t... Is>(std::index_sequence<Is...>) {
            return thrust::make_tuple(&std::get<offset + Is>(tuple)[idx]...);
        };

        // std::make_index_sequence<std::tuple_size_v<decltype(tuple)>>()

        return thrust::make_tuple(
            &v.weight_[idx], &v.charge_[idx],                                  //
            seq.template operator()<2>(std::make_index_sequence<dim>()),       //
            seq.template operator()<2 + dim>(std::make_index_sequence<dim>()), //
            seq.template operator()<2 + dim + dim>(std::make_index_sequence<3>()));

        // return std::tuple_cat(
        //     std::forward_as_tuple(v.weight_, v.charge_),
        //     for_N<dim, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return v.iCell_[i];
        //     }), for_N<dim, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return
        //     v.delta_[i]; }), for_N<3, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return
        //     v.v_[i]; }) //
        // );
    }

    template<typename Particles>
    auto static make(Particles& particles, std::size_t const& idx = 0) _PHARE_ALL_FN_
    {
        auto static constexpr dim = Particles::dimension;

        if constexpr (Particles::storage_mode == StorageMode::SPAN)
        {
            // thrust::tuple_cat

            return thrust::make_zip_iterator(soaxv_tuple_as_thrust_tuple<dim>(particles, idx));
        }

        // else
        // {
        //     return thrust::make_zip_iterator(       //
        //         thrust::make_tuple(                 //
        //             particles.charge_.data() + idx, // 0
        //             particles.weight_.data() + idx, // 1
        //             particles.iCell_.data() + idx,  // 2
        //             particles.delta_.data() + idx,  // 3
        //             particles.v_.data() + idx       // 4
        //             ));
        // }
    }
    template<typename T>
    static auto& charge(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<0>(it);
    }
    template<typename T>
    static auto& weight(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<1>(it);
    }

    template<typename T>
    static auto& iCell(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<2>(it);
    }
    template<typename T>
    static auto& delta(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<3>(it);
    }
    template<typename T>
    static auto& v(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<4>(it);
    }
};


template<typename SoAVXParticles_t>
struct SoAVXZipParticle
{
    auto constexpr static dim = SoAVXParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAVXIteratorAdaptor::template make<SoAVXParticles_t>(
        std::declval<SoAVXParticles_t&>()))>;

    SoAVXZipParticle(SoAVXParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAVXIteratorAdaptor::make(ps, i)}
    {
    }

    auto& charge() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& weight() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& iCell() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    auto& delta() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    auto& v() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }

    Iterator it;
};



template<typename SoAVXParticles_t>
struct SoAVXZipConstParticle
{
    auto constexpr static dim = SoAVXParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAVXIteratorAdaptor::template make<SoAVXParticles_t>(
        std::declval<SoAVXParticles_t&>()))>;

    SoAVXZipConstParticle(SoAVXParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAVXIteratorAdaptor::make(ps, i)}
    {
    }

    auto& charge() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }

    Iterator it;
};


} // namespace PHARE::core::detail


namespace PHARE::core
{

template<typename SoAParticles_t, bool _is_const = false>
struct SoAZipParticle_t
{
    bool static constexpr is_const
        = _is_const || std::is_const_v<std::remove_reference_t<SoAParticles_t>>;

    using value_type = std::conditional_t<is_const, detail::SoAZipConstParticle<SoAParticles_t>,
                                          detail::SoAZipParticle<SoAParticles_t>>;
};

template<typename SoAParticles_t, bool _is_const = false>
struct SoAVXZipParticle_t
{
    bool static constexpr is_const
        = _is_const || std::is_const_v<std::remove_reference_t<SoAParticles_t>>;

    using value_type = std::conditional_t<is_const, detail::SoAVXZipConstParticle<SoAParticles_t>,
                                          detail::SoAVXZipParticle<SoAParticles_t>>;
};


template<typename Particles>
auto particle_zip_iterator(Particles& ps, std::size_t const i)
{
    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, SoA))
        return typename SoAZipParticle_t<Particles>::value_type{ps, i};
    else if constexpr (any_in(Particles::layout_mode, SoAVX))
        return typename SoAVXZipParticle_t<Particles>::value_type{ps, i};
    else
        throw std::runtime_error("no impl");
}

template<typename T, std::size_t dim>
auto partitionner(detail::SoAIteratorAdaptor& begin, detail::SoAIteratorAdaptor& end,
                  Box<T, dim> const& box)
{
    return partitionner(begin, end, box, [&box](auto const& part) {
        return isIn(detail::SoAIteratorAdaptor::iCell(part), box);
    });
}

} // namespace PHARE::core


#endif // PHARE_HAVE_THRUST


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP */
