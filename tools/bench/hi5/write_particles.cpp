
#ifndef PHARE_DIAG_DOUBLES
#define PHARE_DIAG_DOUBLES 0
#endif

#include "diagnostic/detail/h5writer.hpp"
#include "diagnostic/diagnostic_manager.hpp"
#include "hdf5/detail/hdf5_utils.hpp"
#include "hdf5/detail/h5/h5_file.hpp"

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "phare/phare.hpp"

#include "benchmark/benchmark.h"

std::size_t constexpr dim = 1;

namespace PHARE::diagnostic
{

void do_bench(benchmark::State& state)
{
    using HiFile            = HighFive::File;
    using Packer            = core::ParticlePacker<dim>;
    using ParticleArray     = core::ParticleArray<dim>;
    using ParticleArray_SOA = core::ParticleArray<dim, true>;

    auto getSize = [](auto const& value) -> std::size_t {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (hdf5::is_array_dataset<ValueType, dim>)
            return value.size();
        else
            return 1u; /* not an array so value one of type ValueType*/
    };

    auto createDataSet_ = [&](auto& hi5, auto const& path, auto const size, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (hdf5::is_array_dataset<ValueType, dim>)
            return hi5.template create_data_set<typename ValueType::value_type>(
                path, HighFive::DataSpace(size));
        else
            return hi5.template create_data_set<ValueType>(path, HighFive::DataSpace(size));
    };


    auto writeParticles = [](auto& datasets, auto const& particles) {
        datasets[0].write(particles.weight.data());
        datasets[1].write(particles.charge.data());
        datasets[2].write(particles.iCell.data());
        datasets[3].write(particles.delta.data());
        datasets[4].write(particles.v.data());
    };

    auto keys = core::packer_keys();
    std::string path{"/lol/"};
    while (state.KeepRunning())
    {
        hdf5::h5::HighFiveFile hi5("lol.lol",
                                   HiFile::ReadWrite | HiFile::Create | HiFile::Truncate);
        auto d = hi5.create_data_set<float>("/No", 1);
        std::vector<decltype(d)> datasets;

        ParticleArray_SOA particles{100000};
        ParticleArray particleArray(100000);
        Packer{particleArray}.pack(particles);

        std::size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            datasets.emplace_back(
                createDataSet_(hi5, path + keys[part_idx], getSize(arg) * particles.size(), arg));
            part_idx++;
        });
        writeParticles(datasets, particles);
    }
}


} // namespace PHARE::diagnostic


BENCHMARK(PHARE::diagnostic::do_bench)->Unit(benchmark::kMicrosecond);


int main(int argc, char* argv[])
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    return 0;
}
