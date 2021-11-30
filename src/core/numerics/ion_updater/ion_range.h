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
template<typename HybridStateView>
struct SolverPPCUpdateDAO
{
    using ParticleArray   = typename HybridStateView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, HybridStateView>;

    SolverPPCUpdateDAO(ParticleRange_t domain_)
        : domain{domain_}
    {
    }

    ParticleRange_t domain;
    std::unique_ptr<ParticleRange_t> patch_ghost, level_ghost;
};


template<typename HybridStateView, typename Map_t> //
auto pin_ghosts_to_first_range(std::vector<std::vector<Map_t>>& ranges_map_vec)
{
    using ParticleArray   = typename HybridStateView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, HybridStateView>;

    auto make_ghost_range = [](auto& ghosts, auto& dao) {
        return std::make_unique<ParticleRange_t>(makeRange(ghosts), dao.domain.pop_idx,
                                                 dao.domain.view);
    };

    std::unordered_set<ParticleArray*> assigned_patch_ghost;

    for (auto& ranges_maps : ranges_map_vec)
        for (auto& map : ranges_maps)
            for (auto& [_, daos] : map)
            {
                auto& dao = daos[0];
                auto& pop = dao.domain.view->ions[dao.domain.pop_idx];

                if (assigned_patch_ghost.count(pop.patch_ghost) == 0)
                {
                    dao.patch_ghost = make_ghost_range(*pop.patch_ghost, dao);
                    dao.level_ghost = make_ghost_range(*pop.level_ghost, dao);
                    assigned_patch_ghost.emplace(pop.patch_ghost);
                }
            }
}


template<typename HybridStateView, typename ParticleRangesMaps>
auto merge_contiguous_ranges(ParticleRangesMaps& ranges_maps)
{
}


template<typename HybridStateView>
auto solver_update_dao_per_thread(std::vector<std::shared_ptr<HybridStateView>>& views,
                                  std::size_t n_threads)
{
    using ParticleArray   = typename HybridStateView::ParticleArray_t;
    using ParticleRange_t = ParticleRange<typename ParticleArray::iterator, HybridStateView>;
    using Map_t           = std::unordered_map<typename HybridStateView::GridLayout,
                                     std::vector<SolverPPCUpdateDAO<HybridStateView>>>;

    // threads // pops // map<GridLayout, vector<ranges>>
    std::vector<std::vector<Map_t>> maps(n_threads);
    if (views.size() == 0)
        return maps;

    auto n_pops = views[0]->ions.size();

    // pops // threads // ranges
    std::vector<std::vector<std::vector<ParticleRange_t>>> thread_ranges_per_pop
        = PHARE::core::generate(
            [&](auto i) {
                auto particles = PHARE::core::generate(
                    [&](auto& view) { return std::make_tuple(view, view->ions[i].domain); }, views);

                auto accesor = [&](auto& tuple) -> auto& { return *std::get<1>(tuple); };
                auto builder = [&](auto&& range, auto& tuple) {
                    return ParticleRange_t{range, i, std::get<0>(tuple)};
                };

                return PHARE::core::make_balanced_ranges(particles, n_threads, accesor, builder);
            },
            n_pops);

    for (auto& map_vec : maps)
        map_vec.resize(n_pops);

    for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            for (auto& range : thread_ranges_per_pop[pop_idx][t_idx])
                maps[t_idx][pop_idx][range.view->layout].emplace_back(range);

    merge_contiguous_ranges<HybridStateView>(maps);
    pin_ghosts_to_first_range<HybridStateView>(maps);

    return maps;
}

} // namespace PHARE::core


#endif // ION_UPDATER_H
