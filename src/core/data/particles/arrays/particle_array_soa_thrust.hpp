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


#if __has_include(<thrust/iterator/zip_iterator.h>)

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
    // using RefTuple = std::conditional_t<_const_, SoAParticle_crt<dim>, SoAParticle_rt<dim>>

    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)},
          ref{weight(), charge(), iCell(), delta(), v()}
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

    auto& operator*() _PHARE_ALL_FN_ { return ref; }
    auto& operator*() const _PHARE_ALL_FN_ { return ref; }

    Iterator it;
    SoAParticle_rt<dim> ref;
};



template<typename SoAParticles_t>
struct SoAZipConstParticle
{
    auto constexpr static dim = SoAParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipConstParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)},
          ref{weight(), charge(), iCell(), delta(), v()}
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

    auto& operator*() _PHARE_ALL_FN_ { return ref; }
    auto& operator*() const _PHARE_ALL_FN_ { return ref; }

    Iterator it;
    SoAParticle_crt<dim> ref;
};



} // namespace PHARE::core::detail


namespace PHARE::core
{

template<typename T, std::size_t dim>
auto partitionner(detail::SoAIteratorAdaptor& begin, detail::SoAIteratorAdaptor& end,
                  Box<T, dim> const& box)
{
    return partitionner(begin, end, box, [&box](auto const& part) {
        return isIn(detail::SoAIteratorAdaptor::iCell(part), box);
    });
}

} // namespace PHARE::core


#endif // __has_include(<thrust/iterator/zip_iterator.h>)


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP */
