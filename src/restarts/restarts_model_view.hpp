#ifndef RESTART_MODEL_VIEW_HPP
#define RESTART_MODEL_VIEW_HPP

#include "cppdict/include/dict.hpp"


#include "core/utilities/mpi_utils.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "hdf5/phare_hdf5.hpp"


namespace PHARE::restarts
{
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}


template<typename Model>
bool constexpr is_hybrid_model
    = std::is_same_v<solver::type_list_to_hybrid_model_t<typename Model::type_list>, Model>;




template<typename Hierarchy, typename Model, std::enable_if_t<is_hybrid_model<Model>, int> = 0>
class ModelView : public IModelView
{
public:
    ModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }



    auto writeRestartFile(std::string const& path, int timeStepIdx)
    {
        return hierarchy_.writeRestartFile(path, timeStepIdx);
    }


    auto patch_data_ids() { return model_.patch_data_ids(); }


    ModelView(ModelView const&) = delete;
    ModelView(ModelView&&)      = delete;
    ModelView& operator=(ModelView const&) = delete;
    ModelView& operator=(ModelView&&) = delete;

protected:
    Model& model_;
    Hierarchy& hierarchy_;
};



} // namespace PHARE::restarts



#endif // RESTART_MODEL_VIEW_HPP
