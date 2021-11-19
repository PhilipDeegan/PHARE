#ifndef PHARE_BENCH_REAL_BENCH_H
#define PHARE_BENCH_REAL_BENCH_H

#include "mkn/kul/dbg.hpp"
#include "mkn/kul/log.hpp"

#include <atomic>
#include <thread>

#include "core/utilities/types.h"
#include "phare/phare.h"
#include "simulator/simulator.h"
#include "amr/wrappers/hierarchy.h"

#include "bench/core/bench.h"

#include "NumCpp.hpp"

namespace PHARE::real::bench::harris
{
template<typename T, typename SIZE = size_t>
class NCArraySpan : private core::StackVar<nc::NdArray<T>>, public core::Span<T, SIZE>
{
    using Vector = core::StackVar<nc::NdArray<T>>;
    using Span_  = core::Span<T, SIZE>;

public:
    using Vector::var;

    NCArraySpan(nc::NdArray<T>&& that)
        : Vector{std::move(that)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }

    NCArraySpan(std::size_t size, T value = 0)
        : Vector{size, value}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }

    auto& operator()() { return this->var; }
    auto& operator()() const { return this->var; }

    auto static make(std::size_t size, T val = 0)
    {
        auto v = std::make_shared<NCArraySpan<T>>(1, size);
        v->var = val;
        return v;
    }

    auto static make(nc::NdArray<T>&& arr)
    {
        return std::make_shared<NCArraySpan<T>>(std::move(arr));
    }

    auto static asarray(std::vector<T> const& v0)
    {
        return nc::asarray<T>(const_cast<std::vector<T>&>(v0), /*copy=*/false);
    }
};

PHARE::initializer::PHAREDict& register_diags(PHARE::initializer::PHAREDict& dict)
{
    // set before calling this
    // dict["simulation"]["diagnostics"]["filePath"] = std::string{"tools/bench/real"};

    auto& sim          = dict["simulation"];
    auto time_step     = sim["time_step"].template to<double>();
    auto time_step_nbr = sim["time_step_nbr"].template to<int>();
    std::vector<double> timestamps{0, time_step_nbr * time_step};

    std::size_t diag_idx = 0; // per type
    auto& diag_path      = sim["diagnostics"];
    diag_path["mode"]    = std::string{"overwrite"};
    auto add             = [&](std::string type, std::string quantity, std::string prefix = "") {
        auto& type_path                 = diag_path[type];
        auto& name_path                 = type_path[type + std::to_string(diag_idx++)];
        name_path["type"]               = type; //"fluid/electromag/particle";
        name_path["quantity"]           = prefix + quantity;
        name_path["flush_every"]        = std::size_t{1};
        name_path["write_timestamps"]   = timestamps;
        name_path["compute_timestamps"] = timestamps;
        name_path["n_attributes"]       = std::size_t{0};
    };

    diag_idx = 0;
    for (auto qty : {"B", "E"})
        add("electromag", qty, "/EM_");

    // diag_idx = 0;
    // for (int i = 0; i < sim["ions"]["nbrPopulations"].template to<int>(); ++i)
    // {
    //     std::string popName
    //         = sim["ions"]["pop" + std::to_string(i)]["name"].template to<std::string>();
    //     for (auto qty : {"domain"})
    //         add("particle", qty, "/ions/pop/" + popName + "/");
    // }

    // diag_idx = 0;
    // for (int i = 0; i < sim["ions"]["nbrPopulations"].template to<int>(); ++i)
    // {
    //     for (auto qty : {"flux"})
    //         add("fluid", qty, "/ions/pop/" + popName + "/");
    // }
    // for (auto qty : {"density"})
    //     add("fluid", qty, "/ions/");

    return dict;
}


} // namespace PHARE::real::bench::harris

#endif /*PHARE_BENCH_REAL_BENCH_H*/
