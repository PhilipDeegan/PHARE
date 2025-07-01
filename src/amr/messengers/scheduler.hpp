#ifndef PHARE_AMR_MESSENGERS_SCHEDULER_HPP
#define PHARE_AMR_MESSENGERS_SCHEDULER_HPP

#include "refiner.hpp"


#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include <SAMRAI/xfer/VariableFillPattern.h>


namespace PHARE::amr
{


struct Scheduler
{
};


struct RefinerScheduler : public Scheduler
{
    constexpr static std::size_t n_refiner_types = 7; //
    using Schedule                               = SAMRAI::xfer::RefineSchedule;
    using Algorithm                              = SAMRAI::xfer::RefineAlgorithm;



    RefinerScheduler& add(std::unique_ptr<Algorithm> const& algo,
                          std::shared_ptr<Schedule> schedule, int levelNumber)
    {
        schedules_[levelNumber][algo.get()] = schedule;
        return *this;
    }

    template<auto rtype, typename Extra = void>
    RefinerScheduler& add(std::string const& dst,
                          std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                          std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                          std::function<void(double)> optional = {});

    auto& add_algorithm(std::string const dst)
    {
        return algos[dst].emplace_back(std::make_unique<SAMRAI::xfer::RefineAlgorithm>());
    }


    auto& register_resource(auto& rm, auto& dst, auto& src, auto& scratch, auto&&... args)
    {
        auto&& [idDst, idSrc, idScrtch] = rm->getIDsList(dst, src, scratch);
        this->add_algorithm(dst)->registerRefine(idDst, idSrc, idScrtch, args...);
        return *this;
    }


    auto& register_time_interpolated_resource(auto& rm, auto& dst, auto& src, auto& told,
                                              auto& tnew, auto&&... args)
    {
        auto&& [idDst, idSrc, idTold, idTnew] = rm->getIDsList(dst, src, told, tnew);
        this->add_algorithm(dst)->registerRefine(idDst, idSrc, idTold, idTnew, idDst, args...);
        return *this;
    }


    auto& register_vector_field(auto& rm, auto& dst, auto& src, auto& refOp, auto& fillPat)
    {
        return (*this)
            .register_resource(rm, dst.xName, src.xName, dst.xName, refOp, fillPat)
            .register_resource(rm, dst.yName, src.yName, dst.yName, refOp, fillPat)
            .register_resource(rm, dst.zName, src.zName, dst.zName, refOp, fillPat);
    }


    auto& register_time_interpolated_vector_field(auto& rm, auto& dst, auto& src, auto& told,
                                                  auto& tnew, auto&&... args)
    {
        return (*this)
            .register_time_interpolated_resource(rm, dst.xName, src.xName, told.xName, tnew.xName,
                                                 args...)
            .register_time_interpolated_resource(rm, dst.yName, src.yName, told.yName, tnew.yName,
                                                 args...)
            .register_time_interpolated_resource(rm, dst.zName, src.zName, told.zName, tnew.zName,
                                                 args...);
    }


    NO_DISCARD auto& findSchedule(std::unique_ptr<Algorithm> const& algo, int levelNumber) const
    {
        if (!schedules_.count(levelNumber))
            throw std::runtime_error("no schedule for level " + std::to_string(levelNumber));

        if (!schedules_.at(levelNumber).count(algo.get()))
            throw std::runtime_error("Algorithm has not been registered with Communicator");

        return schedules_.at(levelNumber).at(algo.get());
    }


    template<auto rtype>
    void fill(std::string const& dst, int const levelNumber, double const initDataTime)
    {
        assert(bool{phare_funcs[dst][static_cast<int>(rtype)]});
        phare_funcs[dst][static_cast<int>(rtype)](initDataTime);
    }
    template<auto rtype>
    void fill(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime)
    {
        fill<rtype>(dst, level.getLevelNumber(), initDataTime);
    }


    void fill(std::string const& dst, int const levelNumber, double const initDataTime)
    {
        for (auto const& algo : algos[dst])
            this->findSchedule(algo, levelNumber)->fillData(initDataTime);
    }

    void fill(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime)
    {
        fill(dst, level.getLevelNumber(), initDataTime);
    }


    using PhareFuncs
        = std::unordered_map<std::string,
                             std::array<std::function<void(double /*fillTime*/)>, n_refiner_types>>;

    PhareFuncs phare_funcs;
    std::map<int, std::map<Algorithm* const, std::shared_ptr<Schedule>>> schedules_;
    std::unordered_map<std::string, std::vector<std::unique_ptr<Algorithm>>> algos;
};


template<auto Type, typename Extra>
RefinerScheduler& RefinerScheduler::add(
    std::string const& dst, std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
    std::shared_ptr<SAMRAI::hier::PatchLevel> const& level, std::function<void(double)> optional)
{
    auto const levelNumber = level->getLevelNumber();
    if (bool{optional})
        phare_funcs[dst][static_cast<int>(Type)] = optional;
    else
        phare_funcs[dst][static_cast<int>(Type)]
            = [&](double time) { fill(dst, levelNumber, time); };


    for (auto& algo : this->algos[dst])
    {
        // for GhostField we need schedules that take on the level where there is an
        // overlap (there is always for patches lying inside the level) and goes to
        // coarser level where there is not (patch lying on the level border) Create a
        // communication schedule that communicates data within a single level and
        // interpolates data from coarser hierarchy levels where needed.

        // Data will be communicated from the interiors of the source data on the given
        // level to the interiors and ghosts of destination data on the same level where
        // those sources and destinations overlap. Where they do not overlap,
        // data will be interpolated from source data on coarser levels in the patch
        // hierarchy. Data is time interpolated between old and new sources on coarser
        // levels when and where time interpolation is needed and copied from the source
        // components on the patch level into the destination components otherwise. Note
        // that the next coarser level number must correspond to a level in the
        // hierarchy that represents a region of coarser index space than the
        // destination level. Note that the schedule remains valid as long as the levels
        // involved in its creation do not change; thus, it can be used for multiple
        // data communication cycles.
        if constexpr (Type == RefinerType::GhostField)
        {
            this->add(
                algo,
                algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(), hierarchy),
                levelNumber);
        }

        // the following schedule will only fill patch ghost nodes
        // not level border ghosts
        else if constexpr (Type == RefinerType::PatchGhostField)
        {
            this->add(algo, algo->createSchedule(level), levelNumber);
        }

        else if constexpr (Type == RefinerType::PatchFieldBorderSum)
        {
            static_assert(not std::is_same_v<Extra, void> && "NO!");
            this->add(algo,
                      algo->createSchedule(
                          level, 0, std::make_shared<FieldBorderSumTransactionFactory<Extra>>()),
                      levelNumber);
        }

        // this createSchedule overload is used to initialize fields.
        // note that here we must take that createsSchedule() overload and put nullptr
        // as src since we want to take from coarser level everywhere. using the
        // createSchedule overload that takes level, next_coarser_level only would
        // result in interior ghost nodes to be filled with interior of neighbor patches
        // but there is nothing there.
        else if constexpr (Type == RefinerType::InitField)
        {
            this->add(algo, algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                      levelNumber);
        }


        // here we create the schedule that will intialize the particles that lie within
        // the interior of the patches (no ghost, no coarse to fine). We take almost the
        // same overload as for fields above but the version that takes a
        // PatchLevelFillPattern. Here the PatchLevelInteriorFillPattern is used because
        // we want to fill particles only within the interior of the patches of the
        // level. The reason is that filling the their ghost regions with refined
        // particles would not ensure the ghosts to be clones of neighbor patches
        // particles if the splitting from coarser levels is not deterministic.
        else if constexpr (Type == RefinerType::InitInteriorPart)
        {
            this->add(algo,
                      algo->createSchedule(
                          std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(), level,
                          nullptr, levelNumber - 1, hierarchy),
                      levelNumber);
        }

        // here we create a schedule that will refine particles from coarser level and
        // put them into the level coarse to fine boundary. These are the
        // levelGhostParticlesOld particles. we thus take the same createSchedule
        // overload as above but pass it a PatchLevelBorderFillPattern.
        else if constexpr (Type == RefinerType::LevelBorderParticles)
        {
            this->add(
                algo,
                algo->createSchedule(std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(),
                                     level, nullptr, levelNumber - 1, hierarchy),
                levelNumber);
        }


        else if constexpr (Type == RefinerType::ExteriorGhostParticles)
        {
            this->add(algo, algo->createSchedule(level), levelNumber);
        }
    }

    return *this;
}



} // namespace PHARE::amr

#endif
