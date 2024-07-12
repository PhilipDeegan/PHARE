#ifndef PHARE_CORE_DATA_MKN_GPU_HPP
#define PHARE_CORE_DATA_MKN_GPU_HPP

#include "core/def.hpp"

#if PHARE_HAVE_MKN_GPU

#include <vector>
#include "mkn/gpu.hpp"

#include "core/utilities/types.hpp"

namespace PHARE::core::gpu
{

// n threads per cell, 1 block per cell
template<typename Box_t, bool sync = true>
class BoxCellNLauncher : public mkn::gpu::GDLauncher<sync>
{
    using Super                       = mkn::gpu::GDLauncher<sync>;
    static constexpr std::size_t warp = 32; // 32 for nvidia. todo?
public:
    BoxCellNLauncher(Box_t const& box, std::size_t const& _n, size_t const& dev = 0)
        : Super{1, dev}
        , n{_n}
    {
        n_per_cell  = n % warp == 0 ? n : n + (warp - (n % warp));
        this->count = n_per_cell * box.size();
        this->b     = dim3{};
        this->g     = dim3{};
        this->b.x   = n_per_cell;
        this->g.x   = box.size();
    }

    static std::uint32_t block_idx() _PHARE_ALL_FN_ { return blockIdx.x; }
    static std::uint32_t thread_idx() _PHARE_ALL_FN_ { return threadIdx.x; }

private:
    std::size_t n          = 0;
    std::size_t n_per_cell = 0;
};

// one thread per cell
template<typename Box_t, bool sync = true>
class BoxCellLauncher : public mkn::gpu::GDLauncher<sync>
{
    using Super = mkn::gpu::GDLauncher<sync>;

public:
    BoxCellLauncher(Box_t box, size_t dev = 0)
        : Super{box.size(), dev}
    {
    }

private:
};

class AsyncLauncherPool
{
    bool constexpr static sync_on_kernel_finish = false;
    auto static constexpr _accessor             = [](auto& el) -> auto& { return el; };
    // using Accessor                              = decltype(_accessor);

public:
    static auto& INSTANCE()
    {
        static AsyncLauncherPool i;
        return i;
    }

    void check(auto const& containers, std::size_t const& bx = 0)
    {
        launchers.resize(containers.size());
        streams.resize(containers.size());
    }

    template<bool sync = false>
    void operator()(auto&& fn, auto& containers, auto accessor = _accessor)
    {
    }


private:
    std::vector<mkn::gpu::GDLauncher<sync_on_kernel_finish>> launchers;
    std::vector<mkn::gpu::Stream> streams;
};

} // namespace PHARE::core::gpu


#endif // PHARE_HAVE_MKN_GPU
#endif // PHARE_CORE_DATA_MKN_GPU_HPP
