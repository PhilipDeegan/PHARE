#ifndef PHARE_CORE_DATA_FIELD_FIELD_TILES_GPU_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_TILES_GPU_HPP

// device versions of per tile field operations, for GPU_UNIFIED tiled (vec)fields of one patch:
//  one launch each, a block per job, see gpu::FieldTileLauncher
//  each waits for its launch before returning: the jobs live in its scope
//   todo: async (no wait), the jobs must then outlive the launch, see .gpu_todo

#include "core/def.hpp"
#include "core/utilities/kernels.hpp"
#include "core/data/field/field_tiles.hpp"

#include <cstddef>
#include <stdexcept>

namespace PHARE::core
{

namespace detail
{
    // the components of a tiled field or vecfield, for fn(component)
    void for_each_tiled_component(auto& tiled, auto&& fn)
    {
        if constexpr (is_field_tile_set_v<std::remove_cvref_t<decltype(tiled)>>)
            fn(tiled);
        else
            for (auto& c : tiled)
                fn(c);
    }

    template<typename View_t, typename Box_t>
    struct GhostCopyJob
    {
        View_t dst, src;
        Box_t dst_box, src_box;
    };

    // fn(job) on the device, a block per job
    void launch_jobs([[maybe_unused]] auto const& jobs, [[maybe_unused]] auto&& fn)
    {
#if PHARE_HAVE_GPU
        if (jobs.empty())
            return;
        using Launcher = gpu::FieldTileLauncher</*sync=*/true>;
        auto* const jp = jobs.data();
        Launcher{jobs.size()}([=] _PHARE_ALL_FN_() mutable { fn(jp[Launcher::tile_idx()]); });
#else
        throw std::runtime_error("field tiles: GPU_UNIFIED needs a GPU build");
#endif
    }
} // namespace detail


// device version of sync_inner_ghosts: every tile's ghosts from its neighbours' domains, for every
//  component, a block per (tile, neighbour) overlap, threads stride its points
void sync_inner_ghosts_on_device(auto& tiled)
{
#if PHARE_HAVE_GPU
    using Launcher = gpu::FieldTileLauncher</*sync=*/true>;

    auto const jobs = [&]() {
        auto& first      = [&]() -> auto& {
            if constexpr (is_field_tile_set_v<std::remove_cvref_t<decltype(tiled)>>)
                return tiled;
            else
                return tiled[0];
        }();
        using View_t     = std::remove_cvref_t<decltype(first()[0]())>;
        using Box_t      = Box<std::uint32_t, std::remove_cvref_t<decltype(first)>::dimension>;
        using Job_t      = detail::GhostCopyJob<View_t, Box_t>;

        gpu::Vec_t<Job_t> jobs_;
        detail::for_each_tiled_component(tiled, [&](auto& field) {
            for (std::size_t ti = 0; ti < field().size(); ++ti)
                field.for_each_inner_ghost_overlap(
                    ti, [&](auto& t0, auto const& t1, auto const& dst, auto const& src) {
                        jobs_.push_back(Job_t{t0(), t1(), dst, src});
                    });
        });
        return jobs_;
    }();

    using Job_t = std::remove_cvref_t<decltype(jobs)>::value_type;
    detail::launch_jobs(jobs, [] _PHARE_ALL_FN_(Job_t job) {
        Launcher::for_each_index(job.dst_box.size(), [&](auto const i) {
            job.dst(point_at(job.dst_box, i)) = job.src(point_at(job.src_box, i));
        });
    });
#else
    throw std::runtime_error("sync_inner_ghosts_on_device: needs a GPU build");
#endif
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_FIELD_FIELD_TILES_GPU_HPP */
