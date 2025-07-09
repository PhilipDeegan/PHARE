#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "core/logger.hpp"

#if __has_include("bstp/include/BS_thread_pool.hpp")
#define PHARE_HAVE_BSTP_THREAD_POOL 1
#include "bstp/include/BS_thread_pool.hpp"

#else
#error
#define PHARE_HAVE_BSTP_THREAD_POOL 0

#endif

#include <mutex>
#include <memory>
#include <cassert>

namespace PHARE::core
{

struct ThreadPool
{
    ThreadPool()
        : mutices(n_pools)
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        PHARE_LOG_LINE_SS("n_pools " << n_pools);
        PHARE_LOG_LINE_SS("threads_per_pool " << threads_per_pool);
        thread_pools.reserve(n_pools);
        for (std::uint8_t i = 0; i < n_pools; ++i)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
#endif
    }

    static ThreadPool& INSTANCE()
    {
        static ThreadPool i;
        return i;
    }

    auto& get_pool(std::size_t const idx)
    {
        assert(n_pools > 0);
        assert(threads_per_pool > 0);
        assert(idx < thread_pools.size());
        return *thread_pools[idx];
    }

    auto& get_mutex(std::size_t const idx)
    {
        assert(idx < mutices.size());
        return mutices[idx];
    }

    static auto& pool(std::size_t const idx = 0) { return INSTANCE().get_pool(idx); }
    static auto& mutex(std::size_t const idx = 0) { return INSTANCE().get_mutex(idx); }



    // void async(std::size_t const pool_id, auto&& fn)
    // {
    //     assert(thread_pools.size() and pool_id < thread_pools.size());
    //     thread_pools[pool_id]->detach_task(fn);
    // }

    void sync()
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        for (auto& p : thread_pools)
            p->wait();
#endif
    }

    static inline std::size_t n_pools          = 1;
    static inline std::size_t threads_per_pool = 1;

    std::vector<std::mutex> mutices;
#if PHARE_HAVE_BSTP_THREAD_POOL
    std::vector<std::shared_ptr<::BS::thread_pool<::BS::tp::none>>> thread_pools;
#endif
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
