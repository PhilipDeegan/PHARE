
// IS NOT THE SAME AS PYTHON VERSION! :(

#ifndef PHARE_BENCH_REAL_HARRIS_BENCH_H
#define PHARE_BENCH_REAL_HARRIS_BENCH_H

#include "bench/real/real_bench.h"


namespace PHARE::real::bench::harris
{
std::size_t static constexpr dim    = 2;
std::size_t static constexpr interp = 1;
int constexpr cells                 = 10;
auto constexpr dl                   = .2;

using NCArrSpan    = NCArraySpan<double>;
using NCArrSpanPtr = std::shared_ptr<NCArrSpan>;
using Param        = std::vector<double> const&;
using Return       = std::shared_ptr<PHARE::core::Span<double>>;

template<typename Y, typename Y0>
auto S(Y const& y, Y0 const& y0, double const l)
{
    return 0.5 * (1. + nc::tanh((y - y0) / l));
}

template<typename Ret = NCArrSpanPtr>
Ret bx(Param x_in, Param y_in)
{
    auto x  = NCArrSpan::asarray(x_in);
    auto y  = NCArrSpan::asarray(y_in);
    auto Lx = cells, Ly = cells;
    auto w1 = 0.2, w2 = 1.0;
    auto x0 = (x - 0.5 * Lx);
    auto y1 = (y - 0.3 * Ly);
    auto y2 = (y - 0.7 * Ly);
    auto w3 = nc::exp(-(x0 * x0 + y1 * y1) / (w2 * w2));
    auto w4 = nc::exp(-(x0 * x0 + y2 * y2) / (w2 * w2));
    auto w5 = 2.0 * w1 / w2;
    auto v1 = -1., v2 = 1.;

    return NCArrSpan::make(v1 + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
                           + (-w5 * y1 * w3) + (+w5 * y2 * w4));
}

template<typename Ret = NCArrSpanPtr>
Ret by(Param x_in, Param y_in)
{
    auto x = NCArrSpan::asarray(x_in);
    auto y = NCArrSpan::asarray(y_in);

    auto Lx = cells, Ly = cells;
    auto w1 = 0.2, w2 = 1.0;
    auto x0 = (x - 0.5 * Lx);
    auto y1 = (y - 0.3 * Ly);
    auto y2 = (y - 0.7 * Ly);
    auto w3 = nc::exp(-(x0 * x0 + y1 * y1) / (w2 * w2));
    auto w4 = nc::exp(-(x0 * x0 + y2 * y2) / (w2 * w2));
    auto w5 = 2.0 * w1 / w2;

    return NCArrSpan::make((w5 * x0 * w3) + (-w5 * x0 * w4));
}

template<typename Ret = NCArrSpanPtr>
Ret bz(Param x, Param y)
{
    return NCArrSpan::make(x.size());
}

template<typename Ret = NCArrSpanPtr>
Ret b2(Param x, Param y)
{
    auto bx0 = bx(x, y);
    auto by0 = by(x, y);
    auto bz0 = bz(x, y);

    return NCArrSpan::make(nc::power(bx0->var, 2) + nc::power(bx0->var, 2)
                           + nc::power(bz0->var, 2));
}

template<typename Ret = NCArrSpanPtr>
Ret density(Param x0, Param y0)
{
    auto L = cells;
    auto x = NCArrSpan::asarray(x0);
    auto y = NCArrSpan::asarray(y0);

    return NCArrSpan::make(0.2 + 1. / nc::power(nc::cosh((y - L * 0.3) / 0.5), 2)
                           + 1. / nc::power(nc::cosh((y - L * 0.7) / 0.5), 2));
}

auto T(Param x0, Param y0)
{
    auto rho = density(x0, y0);
    auto b2_ = b2(x0, y0);
    auto K   = NCArrSpan::make(x0.size(), 1);

    return NCArrSpan::make(1. / rho->var * (K->var - b2_->var * 0.5));
}

Return vxyz(Param x, Param y)
{
    return NCArrSpan::make(x.size());
}

Return vthxyz(Param x, Param y)
{
    return NCArrSpan::make(nc::sqrt(T(x, y)->var));
}


PHARE::initializer::PHAREDict& createDict()
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;

    auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();

    auto& sim = dict["simulation"];

    sim["dimension"]            = 2;
    sim["interp_order"]         = 1;
    sim["refined_particle_nbr"] = 4;
    sim["time_step"]            = .001;
    sim["time_step_nbr"]        = int{10};
    sim["boundary_type"]        = std::string{"periodic"};

    for (auto s : {"x", "y"})
    {
        sim["grid"]["nbr_cells"][s] = cells;
        sim["grid"]["meshsize"][s]  = dl;
        sim["grid"]["origin"][s]    = 0.;
    }

    auto& amr                              = sim["AMR"];
    amr["max_nbr_levels"]                  = int{2};
    amr["nesting_buffer"]                  = std::vector<int>{0, 0};
    amr["refinement"]["tagging"]["method"] = std::string{"auto"};

    auto& algo                            = sim["algo"];
    algo["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    algo["ohm"]["resistivity"]            = .001;
    algo["ohm"]["hyper_resistivity"]      = .001;

    sim["ions"]["nbrPopulations"] = int{1};
    sim["ions"]["pop0"]["name"]   = std::string{"protons"};
    sim["ions"]["pop0"]["mass"]   = 1.;

    auto& pop_init                 = sim["ions"]["pop0"]["particle_initializer"];
    pop_init["name"]               = std::string{"maxwellian"};
    pop_init["nbr_part_per_cell"]  = int{100};
    pop_init["charge"]             = 1.;
    pop_init["basis"]              = std::string{"cartesian"};
    pop_init["init"]["seed"]       = std::optional<std::size_t>(1337);
    pop_init["density"]            = static_cast<InitFunctionT>(density<Return>);
    pop_init["bulk_velocity_x"]    = static_cast<InitFunctionT>(vxyz);
    pop_init["bulk_velocity_y"]    = static_cast<InitFunctionT>(vxyz);
    pop_init["bulk_velocity_z"]    = static_cast<InitFunctionT>(vxyz);
    pop_init["thermal_velocity_x"] = static_cast<InitFunctionT>(vthxyz);
    pop_init["thermal_velocity_y"] = static_cast<InitFunctionT>(vthxyz);
    pop_init["thermal_velocity_z"] = static_cast<InitFunctionT>(vthxyz);

    sim["electromag"]["name"]             = std::string{"EM"};
    sim["electromag"]["electric"]["name"] = std::string{"E"};
    sim["electromag"]["magnetic"]["name"] = std::string{"B"};

    auto& b_init          = sim["electromag"]["magnetic"]["initializer"];
    b_init["x_component"] = static_cast<InitFunctionT>(bx<Return>);
    b_init["y_component"] = static_cast<InitFunctionT>(by<Return>);
    b_init["z_component"] = static_cast<InitFunctionT>(bz<Return>);

    sim["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    sim["electrons"]["pressure_closure"]["Te"]   = 0.0;

    sim["diagnostics"]["filePath"] = std::string{"tools/bench/real/harris"};
    return register_diags(dict);
}

} // namespace PHARE::real::bench::harris

#endif /*PHARE_BENCH_REAL_HARRIS_BENCH_H*/
