#ifndef PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP


#include "phare_core.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"


#include <cassert>

namespace PHARE::core
{


template<typename ParticleArray_, std::size_t interp_, typename Quantity_ = HybridQuantity>
struct UsableIonsDefaultTypes
{
public:
    auto static constexpr dim         = ParticleArray_::dimension;
    auto static constexpr alloc_mode  = ParticleArray_::alloc_mode;
    auto static constexpr layout_mode = ParticleArray_::layout_mode;
    SimOpts static constexpr opts{dim, interp_, layout_mode, alloc_mode};

    using PHARE_Types     = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t    = PHARE_Types::GridLayout_t;
    using ParticleArray_t = ParticleArray_;
    using Quantity        = Quantity_;

    using UsableVecField_t = UsableVecField<GridLayout_t, alloc_mode, layout_mode>;
    using UsableTensorField_t
        = UsableTensorField<GridLayout_t, /*rank=*/2, alloc_mode, layout_mode>;

    using Grid_t        = UsableTensorField_t::Grid_t;
    using VecField_t    = UsableVecField_t::Super;
    using TensorField_t = UsableTensorField_t::Super;

    using IonPopulation_t = IonPopulation<ParticleArray_t, VecField_t, TensorField_t>;
};


auto inline pop_dict(std::string const& name, std::size_t const ppc = 0)
{
    initializer::PHAREDict popdict;
    popdict["name"]                           = name;
    popdict["mass"]                           = 1.0;
    auto particle_initializer                 = initializer::PHAREDict{};
    particle_initializer["nbr_part_per_cell"] = static_cast<int>(ppc);
    particle_initializer["charge"]            = 1.;
    particle_initializer["basis"]             = std::string{"cartesian"};
    popdict["particle_initializer"]           = particle_initializer;
    return popdict;
}

template<typename _defaults>
class UsableIonsPopulation : public _defaults::IonPopulation_t
{
    using Quantity        = _defaults::Quantity;
    using GridLayout_t    = _defaults::GridLayout_t;
    using VecField_t      = _defaults::VecField_t;
    using TensorField_t   = _defaults::TensorField_t;
    using ParticleArray_t = _defaults::ParticleArray_t;
    using Super           = IonPopulation<ParticleArray_t, VecField_t, TensorField_t>;



    void set()
    {
        auto&& [_F, _M, _pd, _cd, _particles] = Super::getCompileTimeResourcesViewList();
        F.set_on(_F);
        M.set_on(_M);
        _pd.setBuffer(&particleDensity_);
        _cd.setBuffer(&chargeDensity_);
        _particles.setBuffer(&particles.pack());
    }


public:
    UsableIonsPopulation(initializer::PHAREDict const& dict, GridLayout_t const& layout)
        : Super{dict}
        , layout_{layout}
        , particleDensity_{this->name() + "_particleDensity", layout_, Quantity::Scalar::rho, 1}
        , chargeDensity_{this->name() + "_chargeDensity", layout_, Quantity::Scalar::rho, 1}
        , F{this->name() + "_flux", layout, Quantity::Vector::V, 1}
        , M{this->name() + "_momentumTensor", layout, Quantity::Tensor::M, 1}
        , particles{this->name(), layout}
    {
        set();
    }

    UsableIonsPopulation(UsableIonsPopulation const& that)
        : Super{pop_dict(that.name())}
        , layout_{that.layout_}
        , particleDensity_{that.particleDensity_}
        , chargeDensity_{that.chargeDensity_}
        , F{that.F}
        , M{that.M}
        , particles{that.particles}
    {
        set();
    }

    UsableIonsPopulation(UsableIonsPopulation&&)                 = default;
    UsableIonsPopulation& operator=(UsableIonsPopulation const&) = delete;
    UsableIonsPopulation& operator=(UsableIonsPopulation&&)      = default;


    Super& super() { return *this; }
    Super const& super() const { return *this; }
    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    GridLayout_t layout_;
    _defaults::Grid_t particleDensity_, chargeDensity_;
    _defaults::UsableVecField_t F;
    _defaults::UsableTensorField_t M;
    UsableParticlesPopulation<ParticleArray_t> particles;
};


template<typename _defaults>
class UsableIons
    : public Ions<typename _defaults::IonPopulation_t, typename _defaults::GridLayout_t>
{
    using Quantity     = _defaults::Quantity;
    using GridLayout_t = _defaults::GridLayout_t;
    using Super        = Ions<typename _defaults::IonPopulation_t, GridLayout_t>;



    template<typename PopNames>
    auto static super(PopNames const& pop_names, std::size_t const ppc)
    {
        initializer::PHAREDict dict;
        dict["nbrPopulations"] = pop_names.size();
        for (std::size_t i = 0; i < pop_names.size(); ++i)
            dict["pop" + std::to_string(i)] = pop_dict(pop_names[i], ppc);
        return dict;
    }



    auto static super(Super const supe)
    {
        initializer::PHAREDict dict;
        dict["nbrPopulations"] = supe.size();
        for (std::size_t i = 0; i < supe.size(); ++i)
            dict["pop" + std::to_string(i)] = pop_dict(supe[i].name());
        return dict;
    }

    void set()
    {
        auto&& [_bV, _M, _cd, _md] = Super::getCompileTimeResourcesViewList();
        Vi.set_on(_bV);
        M.set_on(_M);
        _cd.setBuffer(&chargeDensity_);
        _md.setBuffer(&massDensity_);

        auto& super_pops = Super::getRunTimeResourcesViewList();
        super_pops.clear();
        for (auto& pop : populations)
            super_pops.emplace_back(*pop);
    }


public:
    using ParticleArray_t = typename _defaults::ParticleArray_t;

    UsableIons(GridLayout_t const& layout, initializer::PHAREDict const& dict)
        : Super{dict}
        , layout_{layout}
        , massDensity_{"massDensity", layout_, Quantity::Scalar::rho, 1}
        , chargeDensity_{"chargeDensity", layout_, Quantity::Scalar::rho, 1}
        , Vi{"bulkVel", layout, Quantity::Vector::V, 1}
        , M{"momentumTensor", layout, Quantity::Tensor::M, 0}
    {
        auto& super_pops = Super::getRunTimeResourcesViewList();
        populations.reserve(super_pops.size());
        for (std::size_t i = 0; i < super_pops.size(); ++i)
            populations.emplace_back(dict["pop" + std::to_string(i)], layout);
        set();
    }

    UsableIons(GridLayout_t const& layout, std::vector<std::string> const& pop_names,
               std::size_t const ppc = 0)
        : UsableIons{layout, super(pop_names, ppc)}
    {
    }

    UsableIons(GridLayout_t const& layout, std::string const& pop_name, std::size_t const ppc = 0)
        : UsableIons{layout, std::vector<std::string>{pop_name}, ppc}
    {
    }

    UsableIons(UsableIons&& that)
        : Super(super(*that))
        , layout_{that.layout_}
        , massDensity_{std::move(that.massDensity_)}
        , chargeDensity_{std::move(that.chargeDensity_)}
        , Vi{std::move(that.Vi)}
        , M{std::move(that.M)}
        , populations{std::move(that.populations)}
    {
        set();
    }

    UsableIons(UsableIons const& that)
        : Super(super(*that))
        , layout_{that.layout_}
        , massDensity_{that.massDensity_}
        , chargeDensity_{that.chargeDensity_}
        , Vi{that.Vi}
        , M{that.M}
        , populations{that.populations}
    {
        set();
    }

    UsableIons& operator=(UsableIons const&) = delete;
    UsableIons& operator=(UsableIons&&)      = default;

    Super& super() { return *this; }
    Super const& super() const { return *this; }
    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    GridLayout_t layout_;
    _defaults::Grid_t massDensity_, chargeDensity_;
    _defaults::UsableVecField_t Vi;
    _defaults::UsableTensorField_t M;
    std::vector<UsableIonsPopulation<_defaults>> populations;
};

template<typename ParticleArray_t, std::size_t interp = 1>
using UsableIonsPopulation_t
    = UsableIonsPopulation<UsableIonsDefaultTypes<ParticleArray_t, interp>>;


template<typename ParticleArray_t, std::size_t interp = 1>
using UsableIons_t = UsableIons<UsableIonsDefaultTypes<ParticleArray_t, interp>>;



} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP */
