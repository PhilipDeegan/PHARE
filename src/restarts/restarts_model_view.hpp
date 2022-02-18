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




template<typename Hierarchy, typename Model>
class ModelView : public IModelView
{
public:
    ModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }



    void writeRestartFile(std::string const& path, int timeStepIdx)
    {
        hierarchy_.writeRestartFile(path, timeStepIdx);
    }



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
