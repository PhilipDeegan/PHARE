#ifndef PHARE_ION_RANGE_H
#define PHARE_ION_RANGE_H

#include <memory>
#include <unordered_set>

#include "core/utilities/range/ranges.h"
#include "core/utilities/range/range_replacer.h"
#include "core/data/particles/particle_range.h"

#include "core/logger.h"


namespace PHARE::core
{
template<typename ParticleRange_t>
struct UpdaterParticleRanges
{
    UpdaterParticleRanges(ParticleRange_t domain_)
        : domain{domain_}
    {
    }

    UpdaterParticleRanges(UpdaterParticleRanges const&) = default;

    ParticleRange_t domain;
    std::shared_ptr<ParticleRange_t> patch_ghost, level_ghost;
};


template<typename ParticleRange_t, typename ThreadRanges> //
auto pin_ghosts_to_first_range(ThreadRanges& thread_ranges_per_pop)
{
    auto make_ghost_range = [](auto& ghosts, auto& dao) {
        return std::make_shared<ParticleRange_t>(makeRange(ghosts), dao.domain.layout,
                                                 dao.domain.em, dao.domain.view);
    };

    std::unordered_set<void*> assigned_patch_ghost;

    // for (auto& ranges_maps : ranges_map_vec)
    //     for (auto& map : ranges_maps)
    //         for (auto& [_, daos] : map)

    for (auto& thread_range : thread_ranges_per_pop)
        for (auto& [layout, em, range] : thread_range)
        {
            auto& dao = daos[0];
            auto& pop = *dao.domain.view;

            if (assigned_patch_ghost.count(pop.patch_ghost) == 0)
            {
                dao.patch_ghost = make_ghost_range(*pop.patch_ghost, dao);
                dao.level_ghost = make_ghost_range(*pop.level_ghost, dao);
                assigned_patch_ghost.emplace(pop.patch_ghost);
            }
        }
}

} // namespace PHARE::core

namespace PHARE::core::detail
{
template<typename ParticleRange, typename GridLayout, typename Electromag, typename IonPopView>
auto make_ranges(
    std::vector<std::tuple<GridLayout, Electromag, std::vector<std::shared_ptr<IonPopView>>>>&
        views_per_patch,
    std::size_t n_threads)
{
    abort_if(views_per_patch.size() == 0);

    // assumes all patches have the same number of populations
    auto n_pops = std::get<2>(views_per_patch[0]).size();

    std::vector<std::vector<std::tuple<GridLayout, Electromag, std::shared_ptr<IonPopView>>>>
        views_per_pop(n_pops);

    for (auto& [layout, em, pops] : views_per_patch)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            views_per_pop[pop_idx].emplace_back(layout, em, pops[pop_idx]);

    return generate(
        [&](auto i) {
            auto particles = generate(
                [&](auto& tuple) {
                    auto& [layout, em, view] = tuple;
                    return std::make_tuple(layout, em, view, view->domain);
                },
                views_per_pop[i]);
            auto accesor = [&](auto& tuple) -> auto& { return *std::get<3>(tuple); };
            auto builder = [&](auto&& range, auto& tuple) {
                auto& [layout, em, view, _] = tuple;
                return ParticleRange{{range}, layout, em, view};
            };
            return PHARE::core::make_balanced_ranges(particles, n_threads, accesor, builder);
        },
        views_per_pop.size());
}

template<typename UpdaterParticleRanges>
auto merge_ranges(std::vector<UpdaterParticleRanges>& ranges)
{
    // using ParticleRange_t = typename UpdaterParticleRanges::ParticleRange_t;

    abort_if(ranges.size() == 0);
    auto contiguous = ranges[0];
    for (std::size_t i = 1; i < ranges.size(); ++i)
    {
        auto& next = ranges[i];
        abort_if(contiguous.domain.end() != next.domain.begin());

        auto new_range = makeRange(contiguous.domain.begin(), next.domain.end());
        contiguous.domain
            = /*ParticleRange_t*/ {new_range, contiguous.domain.pop_idx, contiguous.domain.view};
    }
    return contiguous;
}
}

namespace PHARE::core
{
template<typename GridLayout, typename Electromag, typename IonPopView> // patches // populations
auto updater_ranges_per_thread(
    std::vector<std::tuple<GridLayout, Electromag, std::vector<std::shared_ptr<IonPopView>>>>&
        views,
    std::size_t n_threads = 1)
{
    using ParticleArray = typename IonPopView::particle_array_type;
    using ParticleRange_t
        = ParticleRange<typename ParticleArray::iterator, GridLayout, Electromag, IonPopView>;
    using UpdaterRange_t = UpdaterParticleRanges<ParticleRange_t>;
    using Map_t          = std::unordered_map<GridLayout, std::vector<UpdaterRange_t>>;
    using ThreadView     = std::tuple<GridLayout, Electromag, std::vector<UpdaterRange_t>>;
    using ThreadViews    = std::vector<std::vector<ThreadView>>;

    ThreadViews thread_views(n_threads);
    if (views.size() == 0)
        return thread_views;

    // assumes all patches have the same number of populations
    auto n_pops = std::get<2>(views[0]).size();
    auto maps   = generate([&](auto& i) { return std::vector<Map_t>(n_pops); }, n_threads);
    // auto ions   = generate([&](auto& i) { return GridPairs(n_pops); }, n_threads);
    auto thread_ranges_per_pop = detail::make_ranges<ParticleRange_t>(views, n_threads);

    // for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
    //     for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
    //         for (auto& range : thread_ranges_per_pop[pop_idx][t_idx])
    //             maps[t_idx][pop_idx][range.layout].emplace_back(range);

    pin_ghosts_to_first_range<ParticleRange_t>(thread_ranges_per_pop);

    // for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
    //     for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
    //         for (auto& [grid, daos] : maps[t_idx][pop_idx])
    //             thread_views[t_idx][pop_idx].emplace_back(grid, detail::merge_ranges(daos));

    return thread_ranges_per_pop;
}

} // namespace PHARE::core


#endif // ION_UPDATER_H
