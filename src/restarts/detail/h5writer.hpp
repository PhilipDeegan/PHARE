
#ifndef PHARE_DETAIL_RESTART_HIGHFIVE_HPP
#define PHARE_DETAIL_RESTART_HIGHFIVE_HPP

#include "core/logger.hpp"

#include "restarts/restarts_props.hpp"


namespace PHARE::restarts::h5
{
template<typename ModelView>
class Writer
{
public:
    using This = Writer<ModelView>;

    template<typename Hierarchy, typename Model>
    Writer(Hierarchy& hier, Model& model, std::string const hifivePath)
        : filePath_{hifivePath}
        , modelView_{hier, model}
    {
    }

    ~Writer() {}


    template<typename Hierarchy, typename Model>
    static auto make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
        std::string filePath = dict["filePath"].template to<std::string>();
        return std::make_unique<This>(hier, model, filePath);
    }


    void dump(std::vector<RestartsProperties*> const&, int timestep_idx, double current_timestamp)
    {
        PHARE_LOG_LINE;

        modelView_.writeRestartFile(filePath_, timestep_idx);
    }

    auto& modelView() { return modelView_; }


private:
    double timestamp_ = 0;
    std::string filePath_;
    ModelView modelView_;
};



} // namespace PHARE::restarts::h5

#endif /* PHARE_DETAIL_RESTART_HIGHFIVE_H */
