#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp"

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/CascadePartitioner.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>


#include "initializer/data_provider.hpp"


namespace PHARE::amr
{
template<std::size_t _dimension>
class Integrator
{
    int static constexpr rebalance_coarsest_every_default = 1000;

    bool static _rebalance_coarsest(initializer::PHAREDict const& dict)
    {
        return cppdict::get_value(dict, "simulation/advanced/integrator/rebalance_coarsest", 0) > 0;
    }

    bool static _rebalance_coarsest_on_init(initializer::PHAREDict const& dict)
    {
        return cppdict::get_value(dict, "simulation/advanced/integrator/rebalance_coarsest_on_init",
                                  0)
               > 0;
    }

    std::size_t static _rebalance_coarsest_every(initializer::PHAREDict const& dict)
    {
        auto in
            = cppdict::get_value(dict, "simulation/advanced/integrator/rebalance_coarsest_every",
                                 rebalance_coarsest_every_default);
        if (in < 0)
            throw std::runtime_error("rebalance_coarsest_every must be positive");
        return static_cast<std::size_t>(in);
    }

    bool static _is_tagging_refinement(initializer::PHAREDict const& dict)
    {
        return cppdict::get_value(dict, "simulation/AMR/refinement/tagging/method",
                                  std::string{"none"})
               == std::string{"auto"};
    }

public:
    static constexpr std::size_t dimension = _dimension;

    double advance(double dt)
    {
        bool rebalance_coarsest_now = is_tagging_refinement and rebalance_coarsest
                                      and ((time_step_idx == 0 and rebalance_coarsest_on_init)
                                           or (time_step_idx > 0 and rebalance_coarsest_every > 0
                                               and time_step_idx % rebalance_coarsest_every == 0));

        PHARE_LOG_LINE_STR(is_tagging_refinement
                           << " " << time_step_idx << " " << rebalance_coarsest << " "
                           << rebalance_coarsest_on_init << " " << rebalance_coarsest_every << " "
                           << rebalance_coarsest_now);

        auto new_time = timeRefIntegrator_->advanceHierarchy(dt, rebalance_coarsest_now);
        ++time_step_idx;
        return new_time;
    }

    void initialize() { timeRefIntegrator_->initializeHierarchy(); }

    Integrator(PHARE::initializer::PHAREDict const& dict,
               std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
               std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> timeRefLevelStrategy,
               std::shared_ptr<SAMRAI::mesh::StandardTagAndInitStrategy> tagAndInitStrategy,
               std::shared_ptr<SAMRAI::mesh::CascadePartitioner> loadBalancer, //
               double startTime, double endTime);

private:
    bool is_tagging_refinement                 = false;
    bool rebalance_coarsest                    = false;
    bool rebalance_coarsest_on_init            = false;
    std::size_t time_step_idx                  = 0;
    std::size_t const rebalance_coarsest_every = rebalance_coarsest_every_default;


    std::shared_ptr<SAMRAI::algs::TimeRefinementIntegrator> timeRefIntegrator_;
};




//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------


template<std::size_t dimension>
std::shared_ptr<SAMRAI::tbox::MemoryDatabase>
getUserRefinementBoxesDatabase(PHARE::initializer::PHAREDict const& amr);




template<std::size_t _dimension>
Integrator<_dimension>::Integrator(
    PHARE::initializer::PHAREDict const& dict,
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
    std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> timeRefLevelStrategy,
    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitStrategy> tagAndInitStrategy,
    std::shared_ptr<SAMRAI::mesh::CascadePartitioner> loadBalancer, //
    double startTime, double endTime)
    : is_tagging_refinement{_is_tagging_refinement(dict)}
    , rebalance_coarsest{_rebalance_coarsest(dict)}
    , rebalance_coarsest_on_init{_rebalance_coarsest_on_init(dict)}
    , rebalance_coarsest_every{_rebalance_coarsest_every(dict)}
{
    loadBalancer->setSAMRAI_MPI(
        SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld()); // TODO Is it really needed ?

    auto refineDB    = getUserRefinementBoxesDatabase<dimension>(dict["simulation"]["AMR"]);
    auto standardTag = std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
        "StandardTagAndInitialize", tagAndInitStrategy.get(), refineDB);

    auto clustering = [&]() -> std::shared_ptr<SAMRAI::mesh::BoxGeneratorStrategy> {
        if (!dict["simulation"]["AMR"].contains("clustering"))
            throw std::runtime_error(std::string{"clustering type not specificed"});

        auto clustering_type = dict["simulation"]["AMR"]["clustering"].template to<std::string>();

        if (clustering_type == "berger")
        {
            std::shared_ptr<SAMRAI::tbox::Database> bergerDB
                = std::make_shared<SAMRAI::tbox::MemoryDatabase>("Bergerdb");
            bergerDB->putIntegerVector("max_box_size", std::vector<int>(dimension, 10));
            return std::make_shared<SAMRAI::mesh::BergerRigoutsos>(
                SAMRAI::tbox::Dimension{dimension}, bergerDB);
        }

        if (clustering_type == "tile")
            return std::make_shared<SAMRAI::mesh::TileClustering>(
                SAMRAI::tbox::Dimension{dimension});

        throw std::runtime_error(std::string{"Unknown clustering type "} + clustering_type);
    }();

    auto gridding = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
        hierarchy, "GriddingAlgorithm", std::shared_ptr<SAMRAI::tbox::Database>{}, standardTag,
        clustering, loadBalancer);

    std::shared_ptr<SAMRAI::tbox::Database> db
        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("TRIdb");


    db->putDouble("start_time", startTime);
    db->putDouble("end_time", endTime);
    db->putInteger("max_integrator_steps", 1000000);
    db->putIntegerVector(
        "tag_buffer", std::vector<int>(hierarchy->getMaxNumberOfLevels(),
                                       dict["simulation"]["AMR"]["tag_buffer"].template to<int>()));


    timeRefIntegrator_ = std::make_shared<SAMRAI::algs::TimeRefinementIntegrator>(
        "TimeRefinementIntegrator", db, hierarchy, timeRefLevelStrategy, gridding);
}




template<std::size_t dimension>
std::shared_ptr<SAMRAI::tbox::MemoryDatabase>
getUserRefinementBoxesDatabase(PHARE::initializer::PHAREDict const& amr)
{
    auto const& refinement = amr["refinement"];
    auto maxLevelNumber    = amr["max_nbr_levels"].template to<int>();
    if (refinement.contains("boxes"))
    {
        auto& refDict = refinement["boxes"];

        std::shared_ptr<SAMRAI::tbox::MemoryDatabase> refinementBoxesDatabase
            = std::make_shared<SAMRAI::tbox::MemoryDatabase>("StandardTagAndInitialize");

        // user refinement boxes are always defined at t=0
        // auto at0db = refinementBoxesDatabase->putDatabase("at_0");
        // at0db->putInteger("cycle", 0);
        // auto tag0db = at0db->putDatabase("tag_0");
        std::cout << "tagging method is set to REFINE_BOXES\n";
        refinementBoxesDatabase->putString("tagging_method", "REFINE_BOXES");


        for (int levelNumber = 0; levelNumber < maxLevelNumber; ++levelNumber)
        {
            // not all levels are necessarily specified for refinement
            // cppdict will throw when trying to access key L{i} with i = levelNumber
            std::string levelString{"L" + std::to_string(levelNumber)};
            if (refDict.contains(levelString))
            {
                auto& levelDict = refDict[levelString];
                auto samraiDim  = SAMRAI::tbox::Dimension{dimension};
                auto nbrBoxes   = levelDict["nbr_boxes"].template to<int>();
                auto levelDB
                    = refinementBoxesDatabase->putDatabase("level_" + std::to_string(levelNumber));

                std::vector<SAMRAI::tbox::DatabaseBox> dbBoxes;
                for (int iBox = 0; iBox < nbrBoxes; ++iBox)
                {
                    int lower[dimension];
                    int upper[dimension];
                    auto& boxDict = levelDict["B" + std::to_string(iBox)];

                    lower[0] = boxDict["lower"]["x"].template to<int>();
                    upper[0] = boxDict["upper"]["x"].template to<int>();

                    if constexpr (dimension >= 2)
                    {
                        lower[1] = boxDict["lower"]["y"].template to<int>();
                        upper[1] = boxDict["upper"]["y"].template to<int>();
                    }

                    if constexpr (dimension == 3)
                    {
                        lower[2] = boxDict["lower"]["z"].template to<int>();
                        upper[2] = boxDict["upper"]["z"].template to<int>();
                    }

                    dbBoxes.push_back(SAMRAI::tbox::DatabaseBox(samraiDim, lower, upper));
                }
                levelDB->putDatabaseBoxVector("boxes", dbBoxes);
            }
        } // end loop on levels
        return refinementBoxesDatabase;
    }
    else if (refinement.contains("tagging"))
    {
        std::shared_ptr<SAMRAI::tbox::MemoryDatabase> tagDB
            = std::make_shared<SAMRAI::tbox::MemoryDatabase>("StandardTagAndInitialize");
        tagDB->putString("tagging_method", "GRADIENT_DETECTOR");
        return tagDB;
    }
    return nullptr;
}

} // namespace PHARE::amr


#endif // INTEGRATOR_HPP
