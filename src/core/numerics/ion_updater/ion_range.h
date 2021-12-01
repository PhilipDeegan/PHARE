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
template<typename IonPopView>
struct UpdaterParticleRanges
{
    using ParticleArray   = typename IonPopView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, IonPopView>;

    UpdaterParticleRanges(ParticleRange_t domain_)
        : domain{domain_}
    {
    }

    UpdaterParticleRanges(UpdaterParticleRanges const&) = default;

    ParticleRange_t domain;
    std::shared_ptr<ParticleRange_t> patch_ghost, level_ghost;
};


template<typename IonPopView, typename Map_t> //
auto pin_ghosts_to_first_range(std::vector<std::vector<Map_t>>& ranges_map_vec)
{
    using ParticleArray   = typename IonPopView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, IonPopView>;

    auto make_ghost_range = [](auto& ghosts, auto& dao) {
        return std::make_shared<ParticleRange_t>(makeRange(ghosts), dao.domain.pop_idx,
                                                 dao.domain.view);
    };

    std::unordered_set<ParticleArray*> assigned_patch_ghost;

    for (auto& ranges_maps : ranges_map_vec)
        for (auto& map : ranges_maps)
            for (auto& [_, daos] : map)
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
template<typename ParticleRange, typename IonPopView>
auto make_ranges(std::vector<std::vector<std::shared_ptr<IonPopView>>>& views_per_pop_per_patch,
                 std::size_t n_threads)
{
    abort_if(views_per_pop_per_patch.size() == 0);

    auto n_pops = views_per_pop_per_patch[0]
                      .size(); // assumes all patches have the same number of populations

    std::vector<std::vector<std::shared_ptr<IonPopView>>> views_per_patch_per_pop(n_pops);
    for (std::size_t i = 0; i < views_per_pop_per_patch.size(); ++i)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            views_per_patch_per_pop[pop_idx].emplace_back(views_per_pop_per_patch[i][pop_idx]);

    return generate(
        [&](auto i) {
            auto particles
                = generate([&](auto& view) { return std::make_tuple(view, view->domain); },
                           views_per_patch_per_pop[i]);
            auto accesor = [&](auto& tuple) -> auto& { return *std::get<1>(tuple); };
            auto builder = [&](auto&& range, auto& tuple) {
                return ParticleRange{range, i, std::get<0>(tuple)};
            };
            return PHARE::core::make_balanced_ranges(particles, n_threads, accesor, builder);
        },
        views_per_patch_per_pop.size());
}

template<typename UpdateDAO>
auto merge_ranges(std::vector<UpdateDAO>& ranges)
{
    using ParticleRange_t = typename UpdateDAO::ParticleRange_t;

    abort_if(ranges.size() == 0);
    auto contiguous = ranges[0];
    for (std::size_t i = 1; i < ranges.size(); ++i)
    {
        auto& next = ranges[i];
        abort_if(contiguous.domain.end() != next.domain.begin());

        auto new_range = makeRange(contiguous.domain.begin(), next.domain.end());
        contiguous.domain
            = ParticleRange_t{new_range, contiguous.domain.pop_idx, contiguous.domain.view};
    }
    return contiguous;
}
}

namespace PHARE::core
{
template<typename IonPopView> // patches // populations
auto updater_ranges_per_thread(std::vector<std::vector<std::shared_ptr<IonPopView>>>& views,
                               std::size_t n_threads = 1)
{
    using GridLayout      = typename IonPopView::GridLayout;
    using ParticleArray   = typename IonPopView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, IonPopView>;
    using Map_t    = std::unordered_map<GridLayout, std::vector<UpdaterParticleRanges<IonPopView>>>;
    using GridPair = std::pair<GridLayout, UpdaterParticleRanges<IonPopView>>;
    using GridPairs = std::vector<std::vector<GridPair>>;

    if (views.size() == 0)
        return std::vector<GridPairs>(n_threads);

    auto n_pops = views[0].size(); // assumes all patches have the same number of populations
    auto maps   = generate([&](auto& i) { return std::vector<Map_t>(n_pops); }, n_threads);
    auto ions   = generate([&](auto& i) { return GridPairs(n_pops); }, n_threads);
    auto thread_ranges_per_pop = detail::make_ranges<ParticleRange_t>(views, n_threads);

    for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            for (auto& range : thread_ranges_per_pop[pop_idx][t_idx])
                maps[t_idx][pop_idx][range.view->layout].emplace_back(range);

    pin_ghosts_to_first_range<IonPopView>(maps);

    for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            for (auto& [grid, daos] : maps[t_idx][pop_idx])
                ions[t_idx][pop_idx].emplace_back(grid, detail::merge_ranges(daos));

    return ions;
}

} // namespace PHARE::core


#endif // ION_UPDATER_H
