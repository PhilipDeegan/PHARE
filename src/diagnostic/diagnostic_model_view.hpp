#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"

#include "cppdict/include/dict.hpp"

#include <type_traits>

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}


template<typename Hierarchy, typename Model>
class BaseModelView : public IModelView
{
public:
    using GridLayout = Model::gridlayout_type;
    using VecField   = Model::vecfield_type;
    using PatchProperties
        = cppdict::Dict<float, double, std::size_t, std::vector<int>, std::vector<std::uint32_t>,
                        std::vector<double>, std::vector<std::size_t>, std::string,
                        std::vector<std::string>>;

    BaseModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }

    template<typename Action>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        PHARE::amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                               std::forward<Action>(action), minLevel, maxLevel,
                                               model_);
    }

    NO_DISCARD auto boundaryConditions() const { return hierarchy_.boundaryConditions(); }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return hierarchy_.origin(); }
    NO_DISCARD auto cellWidth() const { return hierarchy_.cellWidth(); }

    NO_DISCARD std::string getLayoutTypeString() const
    {
        return std::string{GridLayout::implT::type};
    }

    NO_DISCARD auto getPatchProperties(std::string patchID, GridLayout const& grid) const
    {
        PatchProperties dict;
        dict["origin"]   = grid.origin().toVector();
        dict["nbrCells"] = core::Point<std::uint32_t, Model::dimension>{grid.nbrCells()}.toVector();
        dict["lower"]    = grid.AMRBox().lower.toVector();
        dict["upper"]    = grid.AMRBox().upper.toVector();
        dict["mpi_rank"] = static_cast<std::size_t>(core::mpi::rank());
        return dict;
    }

    NO_DISCARD static auto getEmptyPatchProperties(PatchProperties dict = {})
    {
        dict["origin"]   = std::vector<double>{};
        dict["nbrCells"] = std::vector<std::uint32_t>{};
        dict["lower"]    = std::vector<int>{};
        dict["upper"]    = std::vector<int>{};
        dict["mpi_rank"] = std::size_t{0};
        return dict;
    }

    NO_DISCARD bool hasTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.count(key);
    }

    NO_DISCARD auto& getTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.at(key);
    }

protected:
    Model& model_;
    Hierarchy& hierarchy_;
};


template<typename Hierarchy, typename Model, typename Enable = void>
class ModelView;


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_hybrid_model_v<Model>>>
    : public BaseModelView<Hierarchy, Model>
{
    using VecField = typename Model::vecfield_type;

public:
    using Model_t = Model;
    using BaseModelView<Hierarchy, Model>::BaseModelView;

    NO_DISCARD std::vector<VecField*> getElectromagFields() const
    {
        return {&this->model_.state.electromag.B, &this->model_.state.electromag.E};
    }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }
};


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_mhd_model_v<Model>>>
    : public BaseModelView<Hierarchy, Model>
{
    using Field    = typename Model::field_type;
    using VecField = typename Model::vecfield_type;

public:
    using Model_t = Model;
    using BaseModelView<Hierarchy, Model>::BaseModelView;
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
