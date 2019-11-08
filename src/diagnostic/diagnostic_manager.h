
#ifndef PHARE_Diagnostic_MANAGER_HPP_
#define PHARE_Diagnostic_MANAGER_HPP_

#include "utilities/types.h"
#include "physical_models/hybrid_model.h"
#include "resources_manager/resources_manager.h"


#include <utility>

namespace PHARE
{
namespace diagnostic
{
    enum class Mode { LIGHT, FULL };

    struct Types
    {
        static std::string_view constexpr rho_s = "rho_s", flux_s = "flux_s", E = "E", B = "B",
                                          space_box = "space_box";
    };

} /*namespace diagnostic*/



struct DiagnosticDAO // DAO = DataAccessObject
{
    size_t compute_every = 1, write_every = 1;
    size_t start_iteration = 0, end_iteration = 100; /* likely to be time rather than index*/
                                                     /* do we allow ranges?*/
    std::string name, species, type;
};



class DiagnosticWriter
{
public:
    virtual void write(DiagnosticDAO&)   = 0;
    virtual void compute(DiagnosticDAO&) = 0;
    virtual ~DiagnosticWriter() {}
};



class IDiagnosticsManager
{
public:
    virtual ~IDiagnosticsManager() {}
    IDiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict<1>& dict);
    IDiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict<1>&& dict)
    {
        return addDiagDict(dict);
    }
    void addDiagnostic(DiagnosticDAO& diagnostic) { diagnostics_.emplace_back(diagnostic); }

    virtual void dump() = 0;

protected:
    std::vector<DiagnosticDAO> diagnostics_;
};




class NoOpDiagnosticManager : public IDiagnosticsManager
{
public:
    void dump() override {}
};




IDiagnosticsManager& IDiagnosticsManager::addDiagDict(PHARE::initializer::PHAREDict<1>& dict)
{
    size_t &compute_every   = dict["compute_every"].template to<std::size_t>(),
           &write_every     = dict["write_every"].template to<std::size_t>(),
           &start_iteration = dict["start_iteration"].template to<std::size_t>(),
           &end_iteration   = dict["end_iteration"].template to<std::size_t>();
    std::string &name       = dict["name"].template to<std::string>(),
                &species    = dict["species"].template to<std::string>(),
                &type       = dict["type"].template to<std::string>();

    diagnostics_.emplace_back(PHARE::core::aggregate_adapter<DiagnosticDAO>(
        compute_every, write_every, start_iteration, end_iteration, name, species, type));
    return *this;
}



template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    DiagnosticsManager(Writer& writer)
        : writer_(writer)
    {
    }

    void dump() override;

private:
    Writer& writer_;

    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};



/*TODO
   iterations
*/
template<typename Writer>
void DiagnosticsManager<Writer>::dump(/*time iteration*/)
{
    size_t iter = 1; // TODO replace with time/iteration idx

    auto active = [](auto const& diag, size_t current, size_t every) {
        return (current >= diag.start_iteration && current <= diag.end_iteration)
               && current % every == 0;
    };

    auto needs_write = [&active](auto const& diag, size_t current) {
        return active(diag, current, diag.write_every);
    };

    auto needs_compute = [&active](auto const& diag, size_t current) {
        return active(diag, current, diag.compute_every);
    };

    std::vector<DiagnosticDAO*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        if (needs_compute(diag, iter))
        {
            writer_.getDiagnosticWriterForType(diag.type)->compute(diag);
        }
        if (needs_write(diag, iter))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }
    writer_.dump(activeDiagnostics);
}



// Generic Template declaration, to override per Concrete model type
template<typename Model, typename ModelParams>
class DiagnosticModelView
{
};



// HybridModel<Args...> specialization
template<typename ModelParams>
class DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>, ModelParams>
{
public:
    using Model      = solver::type_list_to_hybrid_model_t<ModelParams>;
    using VecField   = typename Model::vecfield_type;
    using GridLayout = typename Model::gridLayout_type;
    using Attributes = cppdict::Dict<float, double, size_t, std::string>;

    static constexpr auto dimension = Model::dimension;

    DiagnosticModelView(Model& model)
        : model_{model}
    {
    }

    std::vector<VecField*> getElectromagFields() const
    {
        return {&model_.state.electromag.B, &model_.state.electromag.E};
    }

    auto& getIons() const { return model_.state.ions; }

    auto getPatchAttributes(GridLayout& grid);

protected:
    std::string getPatchOrigin(GridLayout& grid)
    {
        auto& bieldx = model_.state.electromag.B.getComponent(core::Component::X);
        return grid
            .fieldNodeCoordinates(bieldx, grid.origin(),
                                  grid.ghostStartIndex(bieldx, core::Direction::X))
            .str();
    }

protected:
    Model& model_;
};



template<typename ModelParams>
auto DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>,
                         ModelParams>::getPatchAttributes(GridLayout& grid)
{
    Attributes dict;
    // dict["id"]     = std::string("id");
    // dict["float"]  = 0.0f;
    // dict["double"] = 0.0;
    // dict["size_t"] = size_t{0};
    dict["origin"] = getPatchOrigin(grid);
    return dict;
}

template<size_t dim>
class ParticlePacker
{
public:
    ParticlePacker(core::ParticleArray<dim> const& particles)
        : particles_{particles}
        , keys_{{"weight", "charge", "iCell", "delta", "v"}}
    {
    }

    auto get(core::Particle<dim> const& particle) const
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }
    auto get(size_t i) { return get(particles_[i]); }

    auto& keys() const { return keys_; }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

    static auto empty()
    {
        core::Particle<dim> particle;
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

private:
    core::ParticleArray<dim> const& particles_;
    std::array<std::string, 5> keys_;
    size_t it_ = 0;
};



// generic subclass of model specialized superclass
template<typename AMRTypes, typename Model>
class AMRDiagnosticModelView : public DiagnosticModelView<Model, typename Model::type_list>
{
public:
    using Super      = DiagnosticModelView<Model, typename Model::type_list>;
    using ResMan     = typename Model::resources_manager_type;
    using GridLayout = typename Model::gridLayout_type;
    using Guard      = amr::ResourcesGuard<ResMan, Model>;
    using Hierarchy  = typename AMRTypes::hierarchy_t;
    using Patch      = typename AMRTypes::patch_t;
    using Super::model_;
    static constexpr auto dimension = Model::dimension;

    AMRDiagnosticModelView(Hierarchy& hierarchy, Model& model)
        : Super{model}
        , hierarchy_{hierarchy}
    {
    }

    auto guardedGrid(Patch& patch) { return GuardedGrid{patch, Super::model_}; }



    template<typename Action, typename... Args>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                        std::forward<Action>(action), minLevel, maxLevel, model_);
    }



protected:
    struct GuardedGrid
    {
        using Guard = typename AMRDiagnosticModelView<AMRTypes, Model>::Guard;

        GuardedGrid(Patch& patch, Model& model)
            : guard_{model.resourcesManager->setOnPatch(patch, model)}
            , grid_{PHARE::amr::layoutFromPatch<GridLayout>(patch)}
        {
        }

        operator GridLayout&() { return grid_; }

        Guard guard_;
        GridLayout grid_;
    };


private:
    Hierarchy& hierarchy_;

    AMRDiagnosticModelView(const AMRDiagnosticModelView&)             = delete;
    AMRDiagnosticModelView(const AMRDiagnosticModelView&&)            = delete;
    AMRDiagnosticModelView& operator&(const AMRDiagnosticModelView&)  = delete;
    AMRDiagnosticModelView& operator&(const AMRDiagnosticModelView&&) = delete;
};



} // namespace PHARE

#endif //  PHARE_Diagnostic_MANAGER_HPP_
