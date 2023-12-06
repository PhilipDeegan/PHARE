#ifndef PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP
#define PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "amr/solvers/solver.hpp"

#include "core/utilities/thread_pool.hpp"

namespace PHARE::solver
{

template<typename GridLayout>
class FaradayTransformer
{
    using core_type = PHARE::core::Faraday<GridLayout>;

public:
    template<typename ViewStates, typename Accessor>
    void operator()(ViewStates& states, Accessor fn, double dt)
    {
        // not thread safe, would probably need to copy (*this) per thread
        for (auto& state : states)
        {
            auto const& [layout, B, E, Bnew] = fn(state);
            auto _                           = core::SetLayout(&layout, faraday_);
            faraday_(B, E, Bnew, dt);
        }
    }

    core_type faraday_;
};

template<typename GridLayout>
class AmpereTransformer
{
    using core_type = PHARE::core::Ampere<GridLayout>;

public:
};

template<typename GridLayout>
class OhmTransformer
{
    using core_type = PHARE::core::Ohm<GridLayout>;

public:
};

template<typename HybridModel_>
struct HybridPPCModelView : public ISolverModelView
{
    using This             = HybridPPCModelView<HybridModel_>;
    using HybridModel_t    = HybridModel_;
    using HybridModel_args = typename HybridModel_t::type_list::Tuple;
    using IPhysicalModel_t = typename HybridModel_t::Interface;

    using Electrons = std::tuple_element_t<3, HybridModel_args>;

    using level_t         = typename HybridModel_t::amr_types::level_t;
    using Electromag      = typename HybridModel_t::electromag_type;
    using Ions            = typename HybridModel_t::ions_type;
    using ParticleArray_t = typename Ions::particle_array_type;
    using Particle_t      = typename ParticleArray_t::value_type;
    using VecFieldT       = typename HybridModel_t::vecfield_type;
    using GridLayout      = typename HybridModel_t::gridlayout_type;
    using Faraday_t       = FaradayTransformer<GridLayout>;
    using Ampere_t        = AmpereTransformer<GridLayout>;
    using Ohm_t           = OhmTransformer<GridLayout>;

    struct PatchState_t;

    using value_type = PatchState_t;

    template<bool _const_ = false>
    struct iterator;

    HybridPPCModelView(level_t& level, IPhysicalModel_t& model)
        : model_{dynamic_cast<HybridModel_t&>(model)}
    {
        regrid(level, model_);
    }

    void regrid(level_t& level, HybridModel_t& hybridModel);

    auto begin() { return iterator</*const=*/false>{*this}; }
    auto begin() const { return iterator</*const=*/true>{*this}; }

    auto end() { return iterator</*const=*/false>{*this, states.size()}; }
    auto end() const { return iterator</*const=*/true>{*this, states.size()}; }

    auto& model() { return model_; }
    auto& model() const { return model_; }

    HybridModel_t& model_;
    std::vector<core::aggregate_adapter<PatchState_t>> states;
    std::vector<PatchState_t*> state_ptrs;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};
};


template<typename HybridModel>
void HybridPPCModelView<HybridModel>::regrid(level_t& level, HybridModel_t& hybridModel)
{
    auto& hybridState = hybridModel.state;
    auto& rm          = *hybridModel.resourcesManager;

    states.clear();

    for (auto& patch : level)
    {
        {
            auto _ = rm.setOnPatch(*patch, hybridState, electromagPred_, electromagAvg_);
            states.emplace_back(                                 //
                PHARE::amr::layoutFromPatch<GridLayout>(*patch), //
                hybridState.ions,                                //
                hybridState.J,                                   //
                hybridState.electromag,                          //
                electromagPred_,                                 //
                electromagAvg_,                                  //
                hybridState.electrons,                           //
                patch                                            //
            );
        }
        assert(states.back() == true);
        assert(states.back().total_particles > 0);
    }
    state_ptrs.reserve(states.size());
    for (auto& state : states)
        state_ptrs.emplace_back(&state);
    std::sort(state_ptrs.begin(), state_ptrs.end(), [](auto const& a, auto const& b) {
        // assert(a->total_particles > 0);
        // assert(b->total_particles > 0);
        return a->total_particles > b->total_particles;
    });
}


template<typename HybridModel>
template<bool _const_>
struct HybridPPCModelView<HybridModel>::iterator
{
    using View = std::conditional_t<!_const_, HybridPPCModelView<HybridModel>,
                                    HybridPPCModelView<HybridModel> const>;

    iterator(View& view_, std::size_t const& i = 0)
        : view{view_}
        , idx{i}
    {
    }

    auto& operator*()
    {
        assert(view.states[idx] == true);
        return view.states[idx];
    }
    auto& operator*() const
    {
        assert(view.states[idx] == true);
        return view.states[idx];
    }

    auto& operator++()
    {
        idx += 1;
        return *this;
    }

    bool operator!=(iterator const& that) const { return &view != &that.view or idx != that.idx; }

    View& view;
    std::size_t idx = 0;
};

template<typename HybridModel>
struct HybridPPCModelView<HybridModel>::PatchState_t
{
    GridLayout layout;
    Ions ions;
    VecFieldT J;
    Electromag electromag, electromagPred, electromagAvg;
    Electrons electrons;
    std::shared_ptr<SAMRAI::hier::Patch> patch;

    std::size_t const total_particles
        = core::sum_from(ions, [](auto const& pop) { return pop.domainParticles().size(); });

    auto ppc_predictor1_faraday()
    {
        return std::forward_as_tuple(layout, electromag.B, electromag.E, electromagPred.B);
    }
    auto ppc_predictor1_ampere() { return std::forward_as_tuple(layout, electromagPred.B, J); }
    auto ppc_predictor1_ohm()
    {
        auto& Ve = electrons.velocity();
        auto& Ne = electrons.density();
        auto& Pe = electrons.pressure();
        return std::forward_as_tuple(layout, Ne, Ve, Pe, electromagPred.B, J, electromagPred.E);
    }


    auto ppc_predictor2_faraday()
    {
        return std::forward_as_tuple(layout, electromag.B, electromagAvg.E, electromagPred.B);
    }
    auto ppc_predictor2_ampere() { return ppc_predictor1_ampere(); }
    auto ppc_predictor2_ohm() { return ppc_predictor1_ohm(); }


    auto ppc_corrector_faraday()
    {
        return std::forward_as_tuple(layout, electromag.B, electromagAvg.E, electromag.B);
    }
    auto ppc_corrector_ampere() { return std::forward_as_tuple(layout, electromag.B, J); }
    auto ppc_corrector_ohm()
    {
        auto& Ve = electrons.velocity();
        auto& Ne = electrons.density();
        auto& Pe = electrons.pressure();
        return std::forward_as_tuple(layout, Ne, Ve, Pe, electromag.B, J, electromag.E);
    }


    auto ppc_average_B()
    {
        return std::forward_as_tuple(electromag.B, electromagPred.B, electromagAvg.B);
    }
    auto ppc_average_E()
    {
        return std::forward_as_tuple(electromag.E, electromagPred.E, electromagAvg.E);
    }


    auto ppc_update_populations() { return std::forward_as_tuple(layout, ions, electromagAvg); }


    operator bool() const
    {
        using Tuple               = decltype((*this)());
        auto constexpr tuple_size = std::tuple_size_v<Tuple>;
        Tuple tup                 = (*this)();

        return core::for_N_all<tuple_size>([&](auto i) {
            auto const& [k, v] = std::get<i>(tup);
            auto isUsable      = v->isUsable();
            if (!isUsable)
            {
                PHARE_LOG_LINE_STR("NOT USABLE: " << k);
            }
            return isUsable;
        });
    }

    auto operator()() const
    {
        return core::make_named_tuple(                        //
            std::make_pair("J", &J),                          //
            std::make_pair("ions", &ions),                    //
            std::make_pair("electrons", &electrons),          //
            std::make_pair("electromag", &electromag),        //
            std::make_pair("electromagAvg", &electromagAvg),  //
            std::make_pair("electromagPred", &electromagPred) //
        );
    }
};


template<typename HybridModel>
class DefaultSolverPPCThreadingStrategy
{
    using This         = DefaultSolverPPCThreadingStrategy<HybridModel>;
    using View_t       = HybridPPCModelView<HybridModel>;
    using State_t      = typename View_t::PatchState_t;
    using IonUpdater_t = core::IonUpdater<typename View_t::Ions, typename View_t::Electromag,
                                          typename View_t::GridLayout>;

    auto static verify_input(std::size_t n_threads)
    {
        if (n_threads == 0)
            throw std::runtime_error("Threads cannot be zero!");
        return n_threads;
    }


public:
    DefaultSolverPPCThreadingStrategy(std::size_t n_threads)
        : n_threads{verify_input(n_threads)}
        , pool{n_threads - 1}
        , patch_views_per_thread{n_threads}
    {
    }


    void solve(initializer::PHAREDict const& updaterDict, View_t& view, double dt,
               core::UpdaterMode mode)
    {
        // auto patch_view_ptrs
        //     = core::generate([](auto& patch_view) { return &patch_view; }, ion_patch_views);

        // if (n_threads > 1)
        //     std::sort(patch_view_ptrs.begin(), patch_view_ptrs.end(),
        //               [](auto& a, auto& b) { return a->total_particles > b->total_particles; });

        for (auto& ppv : patch_views_per_thread)
            ppv.clear();

        {
            std::size_t thread_idx = 0;
            for (auto const& state_ptr : view.state_ptrs)
            {
                patch_views_per_thread[thread_idx++].push_back(state_ptr);
                if (thread_idx == n_threads)
                    thread_idx = 0;
            }
        }

        auto thread_fn = [&](std::size_t thread_idx) {
            IonUpdater_t ionUpdater{updaterDict};
            for (auto& state_ptr : patch_views_per_thread[thread_idx])
            {
                ionUpdater.updatePopulations(state_ptr->ions, state_ptr->electromag,
                                             state_ptr->layout, dt, mode);
            }
        };

        for (std::size_t i = 1; i < n_threads; ++i)
            pool.submit([&, i]() { thread_fn(i); });

        thread_fn(0);
        pool.wait_for_tasks();
    }

private:
    std::size_t n_threads = 1;
    thread_pool pool;

    std::vector<std::vector<State_t*>> patch_views_per_thread;
};


} // namespace PHARE::solver




#endif /* PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP */
