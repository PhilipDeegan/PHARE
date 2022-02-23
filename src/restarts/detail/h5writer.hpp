
#ifndef PHARE_DETAIL_RESTART_HIGHFIVE_HPP
#define PHARE_DETAIL_RESTART_HIGHFIVE_HPP

#include "core/logger.hpp"

#include "restarts/restarts_props.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

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

        auto restart_file = modelView_.writeRestartFile(filePath_, timestep_idx);

        // write model patch_data_ids to file with highfive

        PHARE::hdf5::h5::HighFiveFile h5File{restart_file};

        auto patch_ids = modelView_.patch_data_ids();
        h5File.create_data_set<int>("/phare/patch/ids", patch_ids.size());
        h5File.write_data_set("/phare/patch/ids", patch_ids);
    }

    auto& modelView() { return modelView_; }


private:
    double timestamp_ = 0;
    std::string filePath_;
    ModelView modelView_;
};



} // namespace PHARE::restarts::h5

#endif /* PHARE_DETAIL_RESTART_HIGHFIVE_H */
