#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP

#include <vector>
#include <cstddef>
#include <numeric>
#include <iterator>
#include <stdexcept>

#include "core/def/phare_mpi.hpp" // needs to be before samrai includes


#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <SAMRAI/tbox/RestartManager.h>
#include "SAMRAI/hier/Transformation.h"


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/vector.hpp"


#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/partitionner/partitionner.hpp"


#include "hdf5/detail/hdf5_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"


namespace PHARE
{
namespace amr
{


    /** @brief ParticlesData is a concrete SAMRAI::hier::PatchData subclass to store Particle data
     *
     * This class encapsulates particle storage known by the module core, and by being derived
     * from PatchData is compatible with the SAMRAI data management system.
     *
     * A ParticlesData encapsulates **three** different particle arrays:
     *
     * - domainParticles : these particles are those for which iCell is within the physical domain
     * of the patch
     *
     * - patchGhostParticles: these particles are located within the ghost layer around the physical
     * domain of the patch. We call the "ghost layer" the layer of ghostCellWidth just outside the
     * physical domain of the patch, on borders that have neighbors patchs of the same level.
     * All the particles in the ghost layer are exact clones of particles located on a neighbor
     * patch of the same level. The ghost particles are getting here when then exit the neighbor
     * patch, and can enter the patch.
     *
     * - levelGhostParticles: these particles are located in a layer just passed the patch
     * boundaries that also are level boundaries. These particles are getting here when there is a
     * particle refinement from a coarser level
     *
     */
    /**
     * @brief The ParticlesData class
     */
    template<typename ParticleArray_>
    class ParticlesData : public SAMRAI::hier::PatchData
    {
        using This          = ParticlesData<ParticleArray_>;
        using Super         = SAMRAI::hier::PatchData;
        using ParticleArray = ParticleArray_;
        using SamBox        = SAMRAI::hier::Box;

        // using Particle_t          = typename ParticleArray::Particle_t;
        static constexpr auto is_host_mem
            = ParticleArray::alloc_mode == AllocatorMode::CPU
              || ParticleArray::alloc_mode == AllocatorMode::GPU_UNIFIED;
        static constexpr auto dim = ParticleArray::dimension;
        // add one cell surrounding ghost box to map particles exiting the ghost layer
        static constexpr int ghostSafeMapLayer = 1;

        auto constexpr construct_particles()
        {
            if constexpr (ParticleArray::is_mapped)
                return ParticleArray{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)};
            else
                return ParticleArray{};
        }

        template<typename T, std::size_t S, typename A>
        static auto _put_vec_size(std::vector<std::array<T, S>, A> const& vec)
        {
            return vec.size() * S;
        }
        template<typename T, typename A>
        static auto _put_vec_size(std::vector<T, A> const& vec)
        {
            return vec.size();
        }
        template<typename T, std::size_t S, typename A>
        static auto _put_vec_data(std::vector<std::array<T, S>, A> const& vec)
        {
            return vec[0].data();
        }
        template<typename T, typename A>
        static auto _put_vec_data(std::vector<T, A> const& vec)
        {
            return vec.data();
        }

        template<typename DB, typename T, typename A>
        void static putVecArrayToDB(DB const& db, std::vector<T, A> const& vec,
                                    std::string const& key)
        {
            auto size = _put_vec_size(vec);
            auto data = _put_vec_data(vec);

            using V = std::decay_t<decltype(*data)>;

            if constexpr (std::is_same_v<V, double>)
                db->putDoubleArray(key, data, size);
            else if constexpr (std::is_same_v<V, int>)
                db->putIntegerArray(key, data, size);
            else
                throw std::runtime_error("ParticlesData::putVecArrayToDB unhandled type");
        }

        template<typename T, std::size_t S, typename A>
        static auto _get_vec_size(std::size_t size, std::vector<std::array<T, S>, A> const&)
        {
            return size / S;
        }
        template<typename T, typename A>
        static auto _get_vec_size(std::size_t size, std::vector<T, A> const&)
        {
            return size;
        }
        template<typename T, std::size_t S, typename A>
        static auto _get_vec_data(std::vector<std::array<T, S>, A>& vec)
        {
            return vec[0].data();
        }
        template<typename T, typename A>
        static auto _get_vec_data(std::vector<T, A>& vec)
        {
            return vec.data();
        }

        template<typename DB, typename T, typename A>
        void static getVecArrayFromDB(DB const& db, std::vector<T, A>& vec, std::string const& key)
        {
            vec.resize(_get_vec_size(db->getArraySize(key), vec));

            auto size = db->getArraySize(key);
            auto data = _get_vec_data(vec);

            using V = std::decay_t<decltype(*data)>;

            if constexpr (std::is_same_v<V, double>)
                db->getDoubleArray(key, data, size);

            else if constexpr (std::is_same_v<V, int>)
                db->getIntegerArray(key, data, size);
            else
            {
                throw std::runtime_error("ParticlesData::getVecArrayFromDB unhandled type");
            }
        }

    public:
        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
                      std::string const& name)

            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , domainParticles{construct_particles()}
            , patchGhostParticles{construct_particles()}
            , levelGhostParticles{construct_particles()}
            , levelGhostParticlesOld{construct_particles()}
            , levelGhostParticlesNew{construct_particles()}
            , pack{name,
                   &domainParticles,
                   &patchGhostParticles,
                   &levelGhostParticles,
                   &levelGhostParticlesOld,
                   &levelGhostParticlesNew}
            , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
            , name_{name}
        {
        }


        auto& name() const { return name_; }


        ParticlesData()                     = delete;
        ParticlesData(ParticlesData const&) = delete;
        ParticlesData(ParticlesData&&)      = default;


        ParticlesData& operator=(ParticlesData const&) = delete;

        // SAMRAI interface
        void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
        {
            using Packer = core::ParticlePacker<ParticleArray>;

            Super::putToRestart(restart_db);

            auto putParticles = [&](std::string const& name, auto& particles) {
                // SAMRAI errors on writing 0 size arrays
                if (particles.size() == 0)
                    return;

                if constexpr (ParticleArray::is_mapped)
                    particles.sortMapping();

                Packer packer{particles};
                [[maybe_unused]] core::SoAParticleArray<dim> soa_;

                auto& soa = [&]() -> auto& {
                    if constexpr (ParticleArray::layout_mode == core::LayoutMode::AoS
                                  || ParticleArray::layout_mode == core::LayoutMode::AoSMapped)
                    {
                        soa_.resize(particles.size());
                        packer.pack(soa_);
                        return soa_;
                    }
                    else
                    {
                        return particles;
                    }
                }();

                std::size_t part_idx = 0;
                core::apply(soa.as_tuple(), [&](auto const& v) {
                    using Vector = std::decay_t<decltype(v)>;
                    using T      = typename Vector::value_type;

                    auto const put_vec = [&](auto const& vec) {
                        putVecArrayToDB(restart_db, vec, name + "_" + packer.keys()[part_idx++]);
                    };

                    if constexpr (is_host_mem)
                        put_vec(v);
                    else
                    {
                        static std::vector<T> put_host_vec;
                        put_host_vec.resize(v.size());
                        PHARE::Vector<T>::copy(put_host_vec, v);
                        put_vec(put_host_vec);
                    }
                });
            };

            putParticles("domainParticles", domainParticles);
            putParticles("patchGhostParticles", patchGhostParticles);
            putParticles("levelGhostParticles", levelGhostParticles);
            putParticles("levelGhostParticlesNew", levelGhostParticlesNew);
            putParticles("levelGhostParticlesOld", levelGhostParticlesOld);
        };


        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            using Packer = core::ParticlePacker<ParticleArray>;

            Super::getFromRestart(restart_db);

            auto getParticles = [&](std::string const name, auto& particles) {
                std::array<bool, Packer::n_keys> const keys_exist = core::generate_from(
                    [&](auto const& key) { return restart_db->keyExists(name + "_" + key); },
                    Packer::keys());

                bool all  = core::all(keys_exist);
                bool none = core::none(keys_exist);
                if (!(all or none))
                    throw std::runtime_error("ParticlesData::getFromRestart has been given an "
                                             "invalid input file, inconsistent state detected");

                if (none) // can't read what doesn't exist
                    return;

                auto n_particles
                    = restart_db->getArraySize(name + "_" + Packer::arbitrarySingleValueKey());
                core::SoAParticleArray<dim> soa{n_particles};

                {
                    std::size_t part_idx = 0;
                    core::apply(soa.as_tuple(), [&](auto& arg) {
                        using Vector = std::decay_t<decltype(arg)>;
                        using T      = typename Vector::value_type;

                        auto& vec = [](auto& v) -> auto& {
                            if constexpr (is_host_mem)
                                return v;
                            else
                            {
                                static std::vector<T> get_host_vec;
                                get_host_vec.clear();
                                return get_host_vec;
                            }
                        }(arg);

                        getVecArrayFromDB(restart_db, vec, name + "_" + Packer::keys()[part_idx++]);

                        if constexpr (not is_host_mem)
                            PHARE::Vector<T>::copy(arg, vec);
                    });
                }

                assert(particles.size() == 0);
                particles.reserve(n_particles);
                for (std::size_t i = 0; i < n_particles; ++i)
                    particles.emplace_back(soa.copy(i));

                particles.check();
            };

            getParticles("domainParticles", domainParticles);
            getParticles("patchGhostParticles", patchGhostParticles);
            getParticles("levelGhostParticles", levelGhostParticles);
            getParticles("levelGhostParticlesNew", levelGhostParticlesNew);
            getParticles("levelGhostParticlesOld", levelGhostParticlesOld);
        }


        /**
         * @brief copy takes a source PatchData and tries to copy its particles
         * in our particle arrays where the source ghostbox and our ghost box overlap
         *
         *
         * this function just hands the source ghost box, our ghost box and their
         * intersection to a private copy_
         *
         * We follow the procedure suggested by SAMRAI PatchDatas. If by anychance
         * the source PatchData was not a ParticleData like we are, we'd give ourselves
         * to its copy2 function, assuming it can copy its source content into us.
         */
        void copy(SAMRAI::hier::PatchData const& source) override
        {
            PHARE_LOG_SCOPE(3, "ParticlesData::copy");

            TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

            // throws if fails
            auto& pSource = dynamic_cast<ParticlesData const&>(source);

            SamBox const& sourceBox  = pSource.getBox();
            SamBox const& myGhostBox = getGhostBox();
            const SamBox intersectionBox{sourceBox * myGhostBox};

            if (!intersectionBox.empty())
            {
                copy_(intersectionBox, pSource);
            }
        }


        /**
         * @brief our copy2 will be called by a PatchData if a ParticleData was
         * given to be copied into another kind of PatchData. Here we chose that
         * copy2 throws unconditiionnally.
         */
        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination) const override
        {
            throw std::runtime_error("Cannot cast");
        }


        /**
         * @brief copy with an overlap. Does the copy as the other overload but this time
         * the copy must account for the intersection with the boxes within the overlap
         * The copy is done between the source patch data and myself
         */
        void copy(SAMRAI::hier::PatchData const& source,
                  SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticlesData::copy with overlap");

            // casts throw on failure
            auto& pSource  = dynamic_cast<ParticlesData const&>(source);
            auto& pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap);

            SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
            SAMRAI::hier::BoxContainer const& boxList = pOverlap.getDestinationBoxContainer();
            for (auto const& overlapBox : boxList)
            {
                copy_(overlapBox, pSource, transformation);
            }
        }


        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
                   [[maybe_unused]] SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            throw std::runtime_error("Cannot cast");
        }



        bool canEstimateStreamSizeFromBox() const override { return false; }




        std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            return countNumberParticlesIn_(pOverlap) * ParticleArray::size_of_particle();
        }




        /**
         * @brief packStream is the function that takes particles from our particles arrays
         * that lie in the boxes of the given overlap, and pack them to a stream.
         *
         * Streaming particles means that we have to take particles with iCell on a local source
         * index space , communicate them, and load them at destination with iCell on a destination
         * local index space. To do that we need to:
         *
         * 1- translate source iCell to source AMR index space
         * 2- Apply the offset to shift this AMR index on top of the destination cells
         * 3- pack and communicate particles
         * 4- move back iCell from the shifted AMR index space to the local destination index space
         *
         * Note that step 2 could be done upon reception of the pack, we chose to do it before.
         *
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::packStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            if (pOverlap.isOverlapEmpty())
            {
                constexpr std::size_t zero = 0;
                stream << zero;
            }
            else
            {
                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                ParticleArray outBuffer;
                pack_(pOverlap, transformation, outBuffer);
                stream << outBuffer.size();
                stream.growBufferAsNeeded();
                this->pack_(stream, outBuffer);
            }
        }

        template<typename ParticleArray_t>
        void pack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t const& outBuffer) const
        {
            if constexpr (ParticleArray::layout_mode == core::LayoutMode::SoA)
                std::apply(
                    [&](auto const&... container) {
                        ((stream.pack(container.data(), outBuffer.size())), ...);
                    },
                    outBuffer.as_tuple());
            else
                stream.pack(outBuffer.data(), outBuffer.size());
        }

        /**
         * @brief unpackStream is the function that unpacks a stream of particles to our particle
         * arrays.
         *
         * We get a stream and an overlap. The overlap contains boxes where to put particles and
         * transformation from source to destination AMR indexes.
         *
         * By convention chosen in patckStream, packed particles have their iCell in our AMR index
         * space. This means that before putting them into our local arrays, we need to apply
         * AMRToLocal() to get the proper shift to apply to them
         *
         */
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::unpackStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            if (!pOverlap.isOverlapEmpty())
            {
                // unpack particles into a particle array
                std::size_t numberParticles = 0;
                stream >> numberParticles;
                ParticleArray particleArray(numberParticles);
                unpack_(stream, particleArray);

                // ok now our goal is to put the particles we have just unpacked
                // into the particleData and in the proper particleArray : interior or ghost

                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
                {
                    // we loop over all boxes in the overlap
                    // we have to first take the intersection of each of these boxes
                    // with our ghostBox. This is where unpacked particles should go.

                    SAMRAI::hier::BoxContainer const& overlapBoxes
                        = pOverlap.getDestinationBoxContainer();

                    auto myBox      = getBox();
                    auto myGhostBox = getGhostBox();

                    for (auto const& overlapBox : overlapBoxes)
                    {
                        // our goal here is :
                        // 1/ to check if each particle is in the intersect of the overlap boxes
                        // and our ghostBox 2/ if yes, check if these particles should go within the
                        // interior array or ghost array
                        auto const intersect = getGhostBox() * overlapBox;

                        for (auto const& particle : particleArray)
                        {
                            if (isInBox(intersect, particle))
                            {
                                if (isInBox(myBox, particle))
                                {
                                    domainParticles.push_back(particle);
                                }
                                else
                                {
                                    patchGhostParticles.push_back(particle);
                                }
                            }
                        } // end species loop
                    } // end box loop
                } // end no rotation
            } // end overlap not empty
        }

        template<typename ParticleArray_t>
        void unpack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t& specie) const
        {
            if constexpr (ParticleArray::layout_mode == core::LayoutMode::SoA)
                std::apply(
                    [&](auto&... container) {
                        ((stream.unpack(container.data(), specie.size())), ...);
                    },
                    specie.as_tuple());
            else
                stream.unpack(specie.data(), specie.size());
        }

        core::ParticlesPack<ParticleArray>* getPointer() { return &pack; }


        // Core interface
        // these particles arrays are public because core module is free to use
        // them easily
        ParticleArray domainParticles;
        ParticleArray patchGhostParticles;

        ParticleArray levelGhostParticles;

        ParticleArray levelGhostParticlesOld;
        ParticleArray levelGhostParticlesNew;

        core::ParticlesPack<ParticleArray> pack;

        // using Iterators_t = typename ParticleArray::iterator;
        // std::vector<core::BoxRange<core::Box<int, dim>, Iterators_t>> domain_partition_iterators;


    private:
        //! interiorLocalBox_ is the box, in local index space, that goes from the first to the last
        //! cell in our patch physical domain, i.e. "from dual physical start index to dual physical
        //! end index"
        SamBox interiorLocalBox_;
        std::string name_;

        template<typename SamraiBox>
        void copyMappedParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData);

        template<typename SamraiBox>
        void copyPartitionedParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData);


        void copy_(SamBox const& overlapBox, ParticlesData const& sourceData)
        {
            if constexpr (ParticleArray::is_mapped)
                copyMappedParticles(overlapBox, sourceData);
            else
                copyPartitionedParticles(overlapBox, sourceData);
        }

        template<typename SamraiBox>
        void copyMappedParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData,
                                 SAMRAI::hier::Transformation const& transformation);

        template<typename SamraiBox>
        void copyPartitionedParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData,
                                      SAMRAI::hier::Transformation const& transformation);

        void copy_(SamBox const& overlapBox, ParticlesData const& sourceData,
                   SAMRAI::hier::Transformation const& transformation)
        {
            if constexpr (ParticleArray::is_mapped)
                copyMappedParticles(overlapBox, sourceData, transformation);
            else
                copyPartitionedParticles(overlapBox, sourceData, transformation);
        }


        /**
         * @brief countNumberParticlesIn_ counts the number of particles that lie
         * within the boxes of an overlap. This function count both patchGhost and
         * domain particles since both could be streamed and we want an upperbound
         * on the number of bytes that could be streamed.
         */
        std::size_t countNumberParticlesIn_(SAMRAI::pdat::CellOverlap const& overlap) const
        {
            PHARE_LOG_SCOPE(3, "ParticleData::countNumberParticlesIn_");

            if (overlap.isOverlapEmpty())
                return 0;

            std::size_t numberParticles = 0;
            auto const& overlapBoxes    = overlap.getDestinationBoxContainer();

            for (auto const& overlapBox : overlapBoxes)
            {
                // we are given boxes from the overlap
                // we want to know how many of our local particles
                // lie in that overlap. Overlap is given in the destination
                // index space (see overlap documentation)
                // so we need to transform that overlap box into our box index space.
                // Since source index space + offset = destination indexspace
                // we need to apply an inverseTransform to the overlapBox.
                SamBox shiftedOverlapBox{overlapBox};
                SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
                transformation.inverseTransform(shiftedOverlapBox);
                auto shiftedOverlapBox_p = phare_box_from<dim>(shiftedOverlapBox);
                if constexpr (ParticleArray::is_mapped)
                    numberParticles += domainParticles.nbr_particles_in(shiftedOverlapBox_p);
                else
                {
                    throw std::runtime_error("fix");
                }
            }

            return numberParticles;
        }

        template<typename ParticleArray_t>
        void packMappedParticles(SAMRAI::pdat::CellOverlap const& overlap,
                                 SAMRAI::hier::Transformation const& transformation,
                                 ParticleArray_t& outBuffer) const;

        template<typename ParticleArray_t>
        void packPartitionedParticles(SAMRAI::pdat::CellOverlap const& overlap,
                                      SAMRAI::hier::Transformation const& transformation,
                                      ParticleArray_t& outBuffer) const;

        template<typename ParticleArray_t>
        void pack_(SAMRAI::pdat::CellOverlap const& overlap,
                   SAMRAI::hier::Transformation const& transformation,
                   ParticleArray_t& outBuffer) const
        {
            PHARE_LOG_SCOPE(3, "ParticleData::pack_");
            // we want to put particles from our domain and patchghost arrays
            // that fall into the intersection box Note that the overlap boxes
            // are not in the same index space as our particles.  the
            // transformation offset goes from OUR index space to the
            // destination space.  Therefore we need to inverse transform the
            // overlap box into our index space, intersect each of them with
            // our ghost box and put export them with the transformation offset

            if constexpr (ParticleArray::is_mapped)
                packMappedParticles(overlap, transformation, outBuffer);
            else
                packPartitionedParticles(overlap, transformation, outBuffer);
        }

        // localpatch
        template<typename SamraiBox>
        void copyParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData) const;
        template<typename SamraiBox>
        void copyParticles(SamraiBox const& overlapBox, ParticlesData const& sourceData,
                           SAMRAI::hier::Transformation const& transformation) const;

        // non-local patch
        template<typename ParticleArray_t>
        void packParticles(SAMRAI::pdat::CellOverlap const& overlap,
                           SAMRAI::hier::Transformation const& transformation,
                           ParticleArray_t& outBuffer) const;
    };


} // namespace amr
} // namespace PHARE

namespace PHARE::amr
{
template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyParticles(
    SamraiBox const& intersectionBox, ParticlesData<ParticleArray> const& sourceData) const
{
}
template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyParticles(
    SamraiBox const& intersectionBox, ParticlesData<ParticleArray> const& sourceData,
    SAMRAI::hier::Transformation const& transformation) const
{
}

template<typename ParticleArray>
template<typename ParticleArray_t>
void ParticlesData<ParticleArray>::packParticles(SAMRAI::pdat::CellOverlap const& overlap,
                                                 SAMRAI::hier::Transformation const& transformation,
                                                 ParticleArray_t& outBuffer) const
{
}


template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyPartitionedParticles(
    SamraiBox const& intersectionBox, ParticlesData<ParticleArray> const& sourceData)
{
    std::array particlesArrays{&sourceData.domainParticles, &sourceData.patchGhostParticles};

    auto myDomainBox = this->getBox();

    // for each particles in the source ghost and domain particle arrays
    // we check if it is in the intersectionBox
    // if it is, is it in my domain box ?
    //      - if so, let's add it to my domain particle array
    //      - if not, let's add it to my ghost particle array
    for (auto const& sourceParticlesArray : particlesArrays)
    {
        for (auto const& particle : *sourceParticlesArray)
        {
            if (isInBox(intersectionBox, particle))
            {
                if (isInBox(myDomainBox, particle))
                {
                    domainParticles.push_back(particle);
                }
                else
                {
                    patchGhostParticles.push_back(particle);
                }
            }
        }
    }
}

template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyMappedParticles(
    SamraiBox const& overlapBox, ParticlesData<ParticleArray> const& sourceData)
{
    auto myDomainBox         = this->getBox();
    auto& srcDomainParticles = sourceData.domainParticles;

    PHARE_LOG_START(3, "ParticleData::copy_ DomainToDomain");

    // first copy particles that fall into our domain array
    // they can come from the source domain or patch ghost
    auto destBox  = myDomainBox * overlapBox;
    auto new_size = domainParticles.size();

    if (!destBox.empty())
    {
        auto destBox_p = phare_box_from<dim>(destBox);
        new_size += srcDomainParticles.nbr_particles_in(destBox_p);
        if (domainParticles.capacity() < new_size)
            domainParticles.reserve(new_size);
        srcDomainParticles.export_particles(destBox_p, domainParticles);
    }

    PHARE_LOG_START(3, "ParticlesData::copy_ DomainToGhosts");
    // Now copy particles from the source domain that fall into
    // our ghost layer. The ghost layer is the result of removing the domain box
    // from the intersection box.
    SAMRAI::hier::BoxContainer ghostLayerBoxes{};
    ghostLayerBoxes.removeIntersections(overlapBox, myDomainBox);

    new_size = patchGhostParticles.size();
    for (auto& selectionBox : ghostLayerBoxes)
    {
        if (!selectionBox.empty())
        {
            auto selectionBox_p = phare_box_from<dim>(selectionBox);
            new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
        }
    }
    if (patchGhostParticles.capacity() < new_size)
        patchGhostParticles.reserve(new_size);


    for (auto const& selectionBox : ghostLayerBoxes)
    {
        if (!selectionBox.empty())
        {
            auto selectionBox_p = phare_box_from<dim>(selectionBox);
            srcDomainParticles.export_particles(selectionBox_p, patchGhostParticles);
        }
    }
    PHARE_LOG_STOP(3, "ParticlesData::copy_ DomainToGhosts");
}

template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyPartitionedParticles(
    SamraiBox const& intersectionBox, ParticlesData<ParticleArray> const& sourceData,
    SAMRAI::hier::Transformation const& transformation)
{
    std::array<decltype(sourceData.domainParticles) const*, 2> particlesArrays{
        &sourceData.domainParticles, &sourceData.patchGhostParticles};

    auto myDomainBox = this->getBox();

    auto offset = transformation.getOffset();

    for (auto const& sourceParticlesArray : particlesArrays)
    {
        for (auto const& particle : *sourceParticlesArray)
        {
            // the particle is only copied if it is in the intersectionBox
            // but before its iCell must be shifted by the transformation offset

            auto newParticle = particle.copy();
            for (auto iDir = 0u; iDir < newParticle.iCell().size(); ++iDir)
            {
                newParticle.iCell()[iDir] += offset[iDir];
            }

            if (isInBox(intersectionBox, newParticle))
            {
                // now we now the particle is in the intersection
                // we need to know whether it is in the domain part of that
                // intersection. If it is not, then it must be in the ghost part


                if (isInBox(myDomainBox, newParticle))
                {
                    domainParticles.push_back(newParticle);
                }
                else
                {
                    patchGhostParticles.push_back(newParticle);
                }
            }
        }
    }
}

template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyMappedParticles(
    SamraiBox const& overlapBox, ParticlesData<ParticleArray> const& sourceData,
    SAMRAI::hier::Transformation const& transformation)
{
    auto myDomainBox         = this->getBox();
    auto& srcDomainParticles = sourceData.domainParticles;

    PHARE_LOG_START(3, "ParticleData::copy_ (transform)");

    // first copy particles that fall into our domain array
    // they can come from the source domain or patch ghost
    auto destBox  = myDomainBox * overlapBox;
    auto new_size = domainParticles.size();
    auto offset   = transformation.getOffset();
    auto offseter = [&](auto const& particle) {
        // we make a copy because we do not want to
        // shift the original particle...
        auto shiftedParticle{particle};
        for (std::size_t idir = 0; idir < dim; ++idir)
        {
            shiftedParticle.iCell()[idir] += offset[idir];
        }
        return shiftedParticle;
    };

    PHARE_LOG_START(3, "DomainToDomain (transform)");
    if (!destBox.empty())
    {
        // we cannot select particles from the intersectDomain box
        // right away. The reason is that the transformation may have
        // a non-zero offset and particle iCells from the source are in
        // the source index space, not in the destination index space
        // therefore we need to first modify the destination box to
        // be in the source index space
        // this is done by applying the INVERSE transformation
        // since a *transformation* is from source to destination.

        transformation.inverseTransform(destBox);
        auto destBox_p = phare_box_from<dim>(destBox);
        new_size += srcDomainParticles.nbr_particles_in(destBox_p);

        if (domainParticles.capacity() < new_size)
            domainParticles.reserve(new_size);
        srcDomainParticles.export_particles(destBox_p, domainParticles, offseter);
    }
    PHARE_LOG_STOP(3, "DomainToDomain (transform)");



    PHARE_LOG_START(3, "DomainToGhosts (transform)");
    // Now copy particles from the source domain and patchghost that fall into
    // our ghost layer. The ghost layer is the result of removing the domain box
    // from the intersection box.
    SAMRAI::hier::BoxContainer ghostLayerBoxes{};
    ghostLayerBoxes.removeIntersections(overlapBox, myDomainBox);

    new_size = patchGhostParticles.size();
    for (auto& selectionBox : ghostLayerBoxes)
    {
        if (!selectionBox.empty())
        {
            transformation.inverseTransform(selectionBox);
            auto selectionBox_p = phare_box_from<dim>(selectionBox);
            new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
        }
    }
    if (patchGhostParticles.capacity() < new_size)
        patchGhostParticles.reserve(new_size);


    // ghostLayer boxes already have been inverse transformed
    // in previous loop, not to do again...
    for (auto const& selectionBox : ghostLayerBoxes)
    {
        if (!selectionBox.empty())
        {
            auto selectionBox_p = phare_box_from<dim>(selectionBox);
            srcDomainParticles.export_particles(selectionBox_p, patchGhostParticles, offseter);
        }
    }

    PHARE_LOG_STOP(3, "DomainToGhosts (transform)");
    PHARE_LOG_STOP(3, "ParticleData::copy_ (transform)");
}

template<typename ParticleArray>
template<typename ParticleArray_t>
void ParticlesData<ParticleArray>::packPartitionedParticles(
    SAMRAI::pdat::CellOverlap const& pOverlap, SAMRAI::hier::Transformation const& transformation,
    ParticleArray_t& outBuffer) const
{
    if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
    {
        SAMRAI::hier::BoxContainer const& boxContainer = pOverlap.getDestinationBoxContainer();

        auto const& sourceGhostBox = getGhostBox();

        // sourceBox + offset = source on destination
        // we are given boxes in the Overlap in destination
        // index space. And we want to select all particles
        // in the ghost source box that lie in this overlapBox
        // we thus need to first shift the sourceGhostBox to the
        // destination index space so that its cells (partly) overlap the one
        // of the given overlap boxes.
        // Then pack_ will take all particles which iCell, shifted by the
        // transformation offset onto the overlapBox index space,
        // lie in the overlap box.
        SamBox transformedSource{sourceGhostBox};
        transformation.transform(transformedSource);

        std::array<ParticleArray const*, 2> particlesArrays{&domainParticles, &patchGhostParticles};

        for (auto const& overlapBox : boxContainer)
        {
            SamBox intersectionBox{transformedSource * overlapBox};

            for (auto const& sourceParticlesArray : particlesArrays)
            {
                for (auto const& particle : *sourceParticlesArray)
                {
                    auto shiftedParticle = particle.copy();
                    auto offset          = transformation.getOffset();
                    for (auto i = 0u; i < dim; ++i)
                    {
                        shiftedParticle.iCell()[i] += offset[i];
                    }
                    if (isInBox(intersectionBox, shiftedParticle))
                    {
                        outBuffer.push_back(shiftedParticle);
                    }
                }
            }
        }
    }
    else
    {
        throw std::runtime_error("Error - rotations not handled in PHARE");
    }



    // std::array<ParticleArray const*, 2> particlesArrays{&domainParticles, &patchGhostParticles};

    // for (auto const& sourceParticlesArray : particlesArrays)
    // {
    //     for (auto const& particle : *sourceParticlesArray)
    //     {
    //         auto shiftedParticle{particle};
    //         auto offset = transformation.getOffset();
    //         for (auto i = 0u; i < dim; ++i)
    //         {
    //             shiftedParticle.iCell[i] += offset[i];
    //         }
    //         if (isInBox(intersectionBox, shiftedParticle))
    //         {
    //             outBuffer.push_back(shiftedParticle);
    //         }
    //     }
    // }

    // for (auto box : overlap.getDestinationBoxContainer())
    // {
    //     transformation.inverseTransform(box);
    //     auto box_range = iterator_for_overlap(box);


    //     outBuffer.insert(outBuffer.end(),
    //                      domainParticles.vector().begin() + box_range.begin().idx(),
    //                      domainParticles.vector().begin() + box_range.end().idx());



    //     // void insert(typename std::vector<Particle_t>::iterator position, iterator first,
    //     //             iterator last)
    //     // {
    //     //     particles_.insert(position, first.idx(), last.idx());
    //     // }
    // }
}

template<typename ParticleArray>
template<typename ParticleArray_t>
void ParticlesData<ParticleArray>::packMappedParticles(
    SAMRAI::pdat::CellOverlap const& overlap, SAMRAI::hier::Transformation const& transformation,
    ParticleArray_t& outBuffer) const
{
    auto overlapBoxes = overlap.getDestinationBoxContainer();
    auto offset       = transformation.getOffset();
    std::size_t size  = 0;
    auto offseter     = [&](auto const& particle) {
        auto shiftedParticle{particle};
        for (std::size_t idir = 0; idir < dim; ++idir)
        {
            shiftedParticle.iCell()[idir] += offset[idir];
        }
        return shiftedParticle;
    };
    for (auto const& box : overlapBoxes)
    {
        auto toTakeFrom{box};
        transformation.inverseTransform(toTakeFrom);
        auto toTakeFrom_p = phare_box_from<dim>(toTakeFrom);
        size += domainParticles.nbr_particles_in(toTakeFrom_p);
    }
    outBuffer.reserve(size);
    for (auto const& box : overlapBoxes)
    {
        auto toTakeFrom{box};
        transformation.inverseTransform(toTakeFrom);
        auto toTakeFrom_p = phare_box_from<dim>(toTakeFrom);
        domainParticles.export_particles(toTakeFrom_p, outBuffer, offseter);
    }
}

} // namespace PHARE::amr

#endif
