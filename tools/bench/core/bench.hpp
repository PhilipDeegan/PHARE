#ifndef PHARE_BENCH_CORE_BENCH_H
#define PHARE_BENCH_CORE_BENCH_H

#include "phare_core.hpp"
#include "benchmark/benchmark.h"


namespace PHARE::core::bench
{
template<std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim>,
                                 typename PHARE::core::HybridQuantity::Scalar>;
template<std::size_t dim>
using VecField
    = PHARE::core::VecField<PHARE::core::NdArrayVector<dim>, typename PHARE::core::HybridQuantity>;


template<std::size_t dim>
PHARE::core::Particle<dim> particle(int icell = 15)
{
    return {/*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ PHARE::core::ConstArray<int, dim>(icell),
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.00001, .00001, .00001}}};
}

template<std::size_t dim>
auto make_particles(std::size_t n_particles)
{
    return PHARE::core::ParticleArray<dim>{n_particles, particle<dim>()};
}

template<typename Particles, typename Point>
void disperse(Particles& particles, Point lo, Point up, std::optional<int> seed = std::nullopt)
{
    auto gen = [&]() {
        if (!seed.has_value())
        {
            std::random_device rd;
            std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }();
    for (std::size_t i = 0; i < Particles::dimension; i++)
    {
        std::uniform_int_distribution<> distrib(lo[i], up[i]);
        for (auto& particle : particles)
            particle.iCell[i] = distrib(gen);
    }
}
template<typename Particles>
void disperse(Particles& particles, std::size_t lo, std::size_t up,
              std::optional<int> seed = std::nullopt)
{
    auto constexpr static dim = Particles::dimension;

    disperse(particles, core::ConstArray<int, dim>(lo), core::ConstArray<int, dim>(up), seed);
}

template<std::size_t dim, typename Box>
auto make_particles(std::size_t ppc, Box disperse_in, std::optional<int> seed = std::nullopt)
{
    auto particles = make_particles<dim>(ppc * disperse_in.size());
    disperse(particles, disperse_in.lower, disperse_in.upper, seed);
    return particles;
}


template<typename GridLayout, typename Quantity, std::size_t dim = GridLayout::dimension>
Field<dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}


template<typename Fn, typename Tuple, size_t... Is>
constexpr auto make_tuple_from_(Fn& f, Tuple const& t, std::integer_sequence<size_t, Is...> const&)
{
    return std::make_tuple(f(std::get<Is>(t))...);
}


template<typename Fn, typename Tuple>
constexpr auto make_tuple_from(Fn&& f, Tuple const& t)
{
    return make_tuple_from_(f, t, std::make_integer_sequence<size_t, std::tuple_size_v<Tuple>>{});
}


template<typename GridLayout, typename Tuple>
auto EB(GridLayout const& layout, Tuple const& tuple)
{
    return make_tuple_from(
        [&](auto const& pair) {
            return std::apply([&](auto k, auto v) { return field(k, v, layout); }, pair);
        },
        tuple);
}

template<typename GridLayout, std::size_t dim = GridLayout::dimension>
auto EM(GridLayout const& layout)
{
    return std::make_tuple(EB(layout, HybridQuantity::E_items()),
                           EB(layout, HybridQuantity::B_items()));
}


template<typename GridLayout, std::size_t dim = GridLayout::dimension>
auto rho(GridLayout const& layout)
{
    return field("rho", HybridQuantity::Scalar::rho, layout);
}


template<typename GridLayout>
class VField : public VecField<GridLayout::dimension>
{
public:
    using Super = VecField<GridLayout::dimension>;

    VField(std::string name, GridLayout const& layout)
        : Super{name, HybridQuantity::Vector::V}
        , name{name}
        , xyz{field(name + "x", HybridQuantity::Scalar::Vx, layout),
              field(name + "y", HybridQuantity::Scalar::Vy, layout),
              field(name + "z", HybridQuantity::Scalar::Vz, layout)}
    {
        Super::setBuffer(name + "_x", &xyz[0]);
        Super::setBuffer(name + "_y", &xyz[1]);
        Super::setBuffer(name + "_z", &xyz[2]);
    }

    template<typename _VF_>
    void set_on(_VF_& vf)
    {
        vf.setBuffer(name + "_x", &xyz[0]);
        vf.setBuffer(name + "_y", &xyz[1]);
        vf.setBuffer(name + "_z", &xyz[2]);
    }

protected:
    std::string name;
    std::array<Field<GridLayout::dimension>, 3> xyz;
};

template<typename GridLayout>
class Flux : public VField<GridLayout>
{
public:
    using Super = VField<GridLayout>;

    Flux(GridLayout const& layout, std::string name = "F")
        : Super{name, layout}
    {
    }
};

template<typename GridLayout>
class BulkV : public VField<GridLayout>
{
public:
    using Super = VField<GridLayout>;

    BulkV(GridLayout const& layout)
        : Super{"bulkVel", layout}

    {
    }
};


template<typename GridLayout>
class Electromag : public PHARE::core::Electromag<VecField<GridLayout::dimension>>
{
public:
    using Super = PHARE::core::Electromag<VecField<GridLayout::dimension>>;

    Electromag(GridLayout const& layout)
        : Super{"EM"}
        , emFields{EM(layout)}
    {
        auto& [E, B]       = emFields;
        auto& [ex, ey, ez] = E;
        auto& [bx, by, bz] = B;

        Super::B.setBuffer("EM_B_x", &bx);
        Super::B.setBuffer("EM_B_y", &by);
        Super::B.setBuffer("EM_B_z", &bz);
        Super::E.setBuffer("EM_E_x", &ex);
        Super::E.setBuffer("EM_E_y", &ey);
        Super::E.setBuffer("EM_E_z", &ez);
    }

private:
    decltype(EM(*static_cast<GridLayout*>(0))) emFields;
};

template<typename Ions, typename... Args>
auto single_pop_ions_from(Args&&... args)
{
    static_assert(sizeof...(args) == 5, "Expected 5 arguments");

    auto const& [popName, rho, bulkV, flux, pack] = std::forward_as_tuple(args...);

    initializer::PHAREDict dict;
    dict["nbrPopulations"] = 1;
    initializer::PHAREDict popdict;
    popdict["name"]                 = std::string{popName};
    popdict["mass"]                 = 1.0;
    popdict["particle_initializer"] = initializer::PHAREDict{};
    dict[popName]                   = popdict;
    Ions ions{dict};
    ions.setBuffer("rho", &rho);
    auto& pop0 = ions.getRunTimeResourcesUserList()[0];
    bulkV.set_on(std::get<0>(ions.getCompileTimeResourcesUserList()));
    flux.set_on(std::get<0>(pop0.getCompileTimeResourcesUserList()));
    ions.getRunTimeResourcesUserList()[0].setBuffer(popName, pack);
    ions.getRunTimeResourcesUserList()[0].setBuffer(popName + "_rho", &rho);
    ions.getRunTimeResourcesUserList()[0].density(); // throws on failure;

    return ions;
}


} // namespace PHARE::core::bench

#endif /*PHARE_BENCH_CORE_BENCH_H*/
