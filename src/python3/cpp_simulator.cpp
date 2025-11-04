

#include "core/utilities/thread_pool.hpp"

#include "python3/cpp_simulator.hpp"


#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME cpp
#endif

namespace PHARE::py
{

auto static const n_pools   = core::get_env_as("PHARE_THREAD_POOLS", std::size_t{1});
auto static const n_threads = core::get_env_as("PHARE_THREADS_PER_POOL", std::size_t{1});

// static init!
bool static const premain = []() {
    core::ThreadPool::n_pools          = n_pools;
    core::ThreadPool::threads_per_pool = n_threads;
    core::ThreadPool::INSTANCE(); // init!
    return true;
}();

} // namespace PHARE::py


namespace PHARE::pydata
{

PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m, pybind11::mod_gil_not_used())
{
    declare_essential(m);

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    core::apply(core::possibleSimulators(), [&](auto const& simType) { declare_all(m, simType); });

    declarePatchData<std::vector<double>, 1>(m, "PatchDataVectorDouble_1D");
    declarePatchData<std::vector<double>, 2>(m, "PatchDataVectorDouble_2D");
    declarePatchData<std::vector<double>, 3>(m, "PatchDataVectorDouble_3D");
}

} // namespace PHARE::pydata
