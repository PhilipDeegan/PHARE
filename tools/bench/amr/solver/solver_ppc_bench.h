#ifndef PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H
#define PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H

#include "bench/core/bench.h"

namespace PHARE::amr::bench
{
template<typename HybridPatchState>
struct HybridPatchStateDecomposer
{
    std::size_t constexpr static range_size = 1e5;

    using View          = typename HybridPatchState::view_t;
    using Batch         = std::vector<View>;
    using GridLayout    = typename HybridPatchState::GridLayout_t;
    using ParticleArray = typename HybridPatchState::ParticleArray_t;

    using ParticleRange = PHARE::core::bench::ParticleRange<typename ParticleArray::iterator, View>;
    using HybridStates  = std::vector<std::unique_ptr<HybridPatchState>>;

    auto static make_view(HybridPatchState& state, std::size_t r = 0, std::size_t t = 0)
    {
        return View{state};
    }

    auto static make_shared_view(HybridPatchState& state, std::size_t r = 0, std::size_t t = 0)
    {
        return std::make_shared<View>(state);
    }

    auto static constexpr accesor = [](auto& state) -> auto& { return state->domain; };
    auto static constexpr builder = [](auto&& range, auto& state) {
        return ParticleRange{range, make_shared_view(*state)};
    };


    template<typename Type>
    auto static decompose(HybridStates& states, std::size_t threads,
                          std::integral_constant<Type, 1> ic)
    {
        assert(false);
        return std::vector<Batch>{};
    }

    template<typename Type>
    auto static decompose(HybridStates& states, std::size_t threads,
                          std::integral_constant<Type, 2> ic)
    {
        // using States = HybridStates;



        std::vector<Batch> batches(threads);

        std::size_t per_thread = 0, modulo = 0;
        {
            std::size_t rows = 0;
            for (auto& state : states)
                rows += state->layout.nbrCells()[1];
            per_thread = rows / threads;
            modulo     = rows % threads;
        }

        std::size_t state_idx = 0;
        for (auto& batch : batches)
        {
            std::size_t rows_taken   = 0;
            std::int64_t rows_remain = per_thread;
            while (rows_remain > 0)
            {
                std::size_t rows_remain_ull = static_cast<std::size_t>(rows_remain);
                std::size_t rows_avail      = states[state_idx]->layout.nbrCells()[1];
                if (rows_avail >= rows_remain_ull)
                {
                    batch.emplace_back(make_view(*states[state_idx], rows_avail, rows_taken));
                    rows_taken += rows_avail;
                    rows_remain -= rows_avail;
                }
                else
                {
                    // todo
                }
            }
            ++state_idx;
        }

        // use with
        // PHARE::core::make_balanced_ranges(states, 10, accesor, builder);
        return batches;
    }

    template<typename Type>
    auto static decompose(HybridStates& states, std::size_t threads,
                          std::integral_constant<Type, 3> ic)
    {
        std::vector<Batch> batches(threads);
        assert(false);
        return batches;
    }

    auto static make_balanced_ranges(HybridStates& states, std::size_t n_threads)
    {
        return PHARE::core::make_balanced_ranges(states, n_threads, accesor, builder);
    }
}; // namespace PHARE::amr::bench


template<typename HybridPatchState>
auto decompose(std::vector<std::unique_ptr<HybridPatchState>>& states, std::size_t threads)
{
    return HybridPatchStateDecomposer<HybridPatchState>::decompose(
        states, threads, std::integral_constant<std::size_t, HybridPatchState::dimension>{});
}


template<typename HybridPatchState>
struct SolverPPCUpdateDAO
{
    using ParticleArray   = typename HybridPatchState::ParticleArray_t;
    using ParticleRange_t = PHARE::core::bench::ParticleRange<typename ParticleArray::iterator,
                                                              typename HybridPatchState::view_t>;


    SolverPPCUpdateDAO(ParticleRange_t domain_)
        : domain{domain_}
    {
    }



    ParticleRange_t domain;
    std::unique_ptr<ParticleRange_t> patch_ghost, level_ghost;
};


template<typename HybridPatchState, typename ParticleRangesMaps>
auto pin_ghosts_to_max_operator(ParticleRangesMaps& ranges_maps)
{
    using ParticleArray   = typename HybridPatchState::ParticleArray_t;
    using ParticleRange_t = PHARE::core::bench::ParticleRange<typename ParticleArray::iterator,
                                                              typename HybridPatchState::view_t>;
    using Counts
        = std::unordered_map<typename HybridPatchState::IonPopulation_t const*, std::size_t>;

    auto sorted = [](auto const& map) {
        std::vector<typename Counts::value_type const*> pairs(map.size());
        std::size_t pair_idx = 0;
        for (auto const& pair : map)
            pairs[pair_idx++] = &pair;
        std::sort(pairs.begin(), pairs.end(), // biggest first
                  [](auto const& x, auto const& y) { return x->second > y->second; });
        return pairs;
    };

    std::set<ParticleArray*> assigned_patch_ghost, assigned_level_ghost;
    for (auto& map : ranges_maps)
        for (auto& [_, daos] : map)
            for (auto& dao : daos)
            {
                auto check_set = [&](auto& ptr, auto& ghosts, auto& assigned_ghosts) {
                    if (assigned_ghosts.count(&ghosts) == 0)
                    {
                        ptr = std::make_unique<ParticleRange_t>(makeRange(ghosts), dao.domain.view);
                        assigned_ghosts.insert(&ghosts);
                    }
                };

                Counts counts;
                for (auto const& pop : *dao.domain.view->ions)
                    counts[&pop] += dao.domain.size();
                auto pairs = sorted(counts);

                for (auto const*& pair : pairs)
                {
                    auto const& [pop, _] = *pair;
                    check_set(dao.patch_ghost, pop->patchGhostParticles(), assigned_patch_ghost);
                    check_set(dao.level_ghost, pop->levelGhostParticles(), assigned_level_ghost);
                }
            }
}


template<typename HybridPatchState, typename ParticleRangesMaps>
auto merge_contiguous_ranges(ParticleRangesMaps& ranges_maps)
{
}


template<typename HybridPatchState,
         typename Decomposer   = HybridPatchStateDecomposer<HybridPatchState>,
         typename HybridStates = typename Decomposer::HybridStates>
auto solver_update_dao_per_thread(HybridStates& states, std::size_t n_threads)
{
    using Map_t = std::unordered_map<typename Decomposer::GridLayout,
                                     std::vector<SolverPPCUpdateDAO<HybridPatchState>>>;

    auto ranges_list = Decomposer::make_balanced_ranges(states, n_threads);

    std::vector<Map_t> daos(n_threads);
    for (std::size_t i = 0; i < ranges_list.size(); ++i)
        for (auto const& range : ranges_list[i])
            daos[i][range.view->layout].emplace_back(range);

    merge_contiguous_ranges<HybridPatchState>(daos);
    pin_ghosts_to_max_operator<HybridPatchState>(daos);

    return daos;
}


} // namespace PHARE::amr::bench

#endif /*PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H*/
