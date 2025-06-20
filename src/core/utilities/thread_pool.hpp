#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "core/logger.hpp"

#if __has_include("bstp/include/BS_thread_pool.hpp")
#define PHARE_HAVE_BSTP_THREAD_POOL 1
#include "bstp/include/BS_thread_pool.hpp"

#else
#define PHARE_HAVE_BSTP_THREAD_POOL 0

#endif

#include <memory>
#include <cassert>

namespace PHARE::core
{

struct ThreadPool
{
    ThreadPool()
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        PHARE_LOG_LINE_SS(ThreadPool::threads_per_pool);
        thread_pools.reserve(n_pools);
        for (std::uint8_t i = 0; i < n_pools; i++)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
#endif
    }

    static ThreadPool& INSTANCE()
    {
        static ThreadPool i;
        return i;
    }


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

#if PHARE_HAVE_BSTP_THREAD_POOL
    std::vector<std::shared_ptr<::BS::thread_pool<::BS::tp::none>>> thread_pools;
#endif
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
