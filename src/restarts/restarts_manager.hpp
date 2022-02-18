
#ifndef PHARE_RESTART_MANAGER_HPP_
#define PHARE_RESTART_MANAGER_HPP_

#include "core/logger.hpp"
#include "core/data/particles/particle_array.hpp"

#include "initializer/data_provider.hpp"

#include "restarts_props.hpp"


#include <cmath>
#include <memory>
#include <utility>



namespace PHARE::restarts
{
template<typename RestartManager>
void registerRestarts(RestartManager& rMan, initializer::PHAREDict const& params)
{
    rMan.addRestartDict(params);
}




class IRestartsManager
{
public:
    virtual bool dump(int timeStepIdx, double timeStamp, double timeStep) = 0;
    inline virtual ~IRestartsManager();
};
IRestartsManager::~IRestartsManager() {}



template<typename Writer>
class RestartsManager : public IRestartsManager
{
public:
    bool dump(int timeStepIdx, double timeStamp, double timeStep) override;



    RestartsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::move(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: RestartsManager received null Writer");
    }


    template<typename Hierarchy, typename Model>
    static std::unique_ptr<RestartsManager> make_unique(Hierarchy& hier, Model& model,
                                                        initializer::PHAREDict const& dict)
    {
        auto rMan = std::make_unique<RestartsManager>(Writer::make_unique(hier, model, dict));
        registerRestarts(*rMan, dict);
        return rMan;
    }



    RestartsManager& addRestartDict(initializer::PHAREDict const& dict);
    RestartsManager& addRestartDict(initializer::PHAREDict&& dict) { return addRestartDict(dict); }


    auto& restarts() const { return restarts_; }


    Writer& writer() { return *writer_.get(); }


    RestartsManager(RestartsManager const&) = delete;
    RestartsManager(RestartsManager&&)      = delete;
    RestartsManager& operator=(RestartsManager const&) = delete;
    RestartsManager& operator=(RestartsManager&&) = delete;

private:
    bool needsAction_(double nextTime, double timeStamp, double timeStep)
    {
        // casting to float to truncate double to avoid trailing imprecision
        return static_cast<float>(std::abs(nextTime - timeStamp)) < static_cast<float>(timeStep);
    }


    bool needsWrite_(RestartsProperties const& rest, double const timeStamp, double const timeStep)
    {
        auto const& nextWrite = nextWrite_;

        return nextWrite < rest.writeTimestamps.size()
               and needsAction_(rest.writeTimestamps[nextWrite], timeStamp, timeStep);
    }


    std::vector<RestartsProperties> restarts_;
    std::unique_ptr<Writer> writer_;
    std::size_t nextWrite_ = 0;
};



template<typename Writer>
RestartsManager<Writer>&
RestartsManager<Writer>::addRestartDict(initializer::PHAREDict const& params)
{
    PHARE_LOG_LINE_STR(params.contains("write_timestamps"));

    if (params.contains("write_timestamps"))
    {
        auto& properties           = restarts_.emplace_back();
        properties.writeTimestamps = params["write_timestamps"].template to<std::vector<double>>();
        PHARE_LOG_LINE_STR(properties.writeTimestamps.size());
    }

    PHARE_LOG_LINE_STR(restarts_.size());

    return *this;
}




template<typename Writer>
bool RestartsManager<Writer>::dump(int timeStepIdx, double timeStamp, double timeStep)
{
    PHARE_LOG_LINE;
    std::vector<RestartsProperties*> activeRestarts;

    for (auto& rest : restarts_)
        if (needsWrite_(rest, timeStamp, timeStep))
            activeRestarts.emplace_back(&rest);


    if (activeRestarts.size() > 0)
        writer_->dump(activeRestarts, timeStepIdx, timeStamp);


    if (activeRestarts.size() > 1)
        throw std::runtime_error("Restarts invalid");

    for (auto const* rest : activeRestarts)
        ++nextWrite_;


    return activeRestarts.size() > 0;
}

} // namespace PHARE::restarts

#endif /* PHARE_RESTART_MANAGER_HPP_ */
