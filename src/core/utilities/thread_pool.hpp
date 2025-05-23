#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP


#include "core/logger.hpp"
#include "bstp/include/BS_thread_pool.hpp"

#include <memory>
#include <cassert>

namespace PHARE::core
{

struct ThreadPool
{
    ThreadPool()
    {
        PHARE_LOG_LINE_SS(ThreadPool::threads_per_pool);
        thread_pools.reserve(n_pools);
        for (std::uint8_t i = 0; i < n_pools; i++)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
    }

    static ThreadPool& INSTANCE()
    {
        static ThreadPool i;
        return i;
    }


    void async(std::size_t const pool_id, auto&& fn)
    {
        assert(thread_pools.size() and pool_id < thread_pools.size());
        thread_pools[pool_id]->detach_task(fn);
    }

    void sync()
    {
        for (auto& p : thread_pools)
            p->wait();
    }

    static inline std::size_t n_pools          = 1;
    static inline std::size_t threads_per_pool = 1;

    std::vector<std::shared_ptr<::BS::thread_pool<::BS::tp::none>>> thread_pools;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
