#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

// #include "core/def.hpp"
// #include "core/logger.hpp"


#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"
// #include "core/data/particles/particle_array_exporter.hpp"
#include "core/vector.hpp"
#include "core/utilities/types.hpp"

#include "core/utilities/point/point.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"
// #include "core/utilities/partitionner/partitionner.hpp"

// #include "hdf5/detail/hdf5_utils.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include <amr/data/particles/particles_variable_fill_pattern.hpp>


#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/RestartManager.h>
#include "SAMRAI/hier/Transformation.h"
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <tuple>
#include <vector>
#include <cstddef>
#include <stdexcept>


namespace PHARE
{
namespace amr
{


    /** @brief ParticlesData is a concrete SAMRAI::hier::PatchData subclass
     * to store Particle data
     *
     * A ParticlesData encapsulates **three** different particle arrays:
     *
     * - domainParticles : these particles are those for which iCell
     *   is within the physical domain of the patch
     *
     * - patchGhostParticles: represents particles that left the patch domain and are
     *   physically located in the patch ghost layer of a patch.
     *
     * - levelGhostParticles: represent particles obtained from refinement and
     *   located in level ghost layer. These particles are to be pushed and injected
     *   in domain if they arrive in there.
     *
     *- levelGhostParticlesOld: same as levelGhostParticles but defined at previous
     *  next coarse time step. Used to deposit contribution of these particles
     *  to moments in level ghost nodes
     *
     *- levelGhostParticlesNew: same as levelGhostParticles but defined at next
     * coarser future time step. Used to deposit contribution of these particles
     * to moments in level ghost nodes
     *
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



        auto construct_particles(auto const& ghosts) const
        {
            if constexpr (ParticleArray::is_mapped)
                return ParticleArray{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)};

            else
                return make_particles<ParticleArray>(phare_box_from<dim>(getBox()), ghosts);
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

        void validate_ghosts(auto const ghost)
        {
            core::for_N<dim>([&](auto i) {
                if (ghost[i] != ghost[0])
                    throw std::runtime_error("invalid");
            });
        }

    public:
        // template<typename ParticlesInitializer>
        // ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
        //               std::string const& name, ParticlesInitializer const initializer)
        //     : SAMRAI::hier::PatchData::PatchData(box, ghost)
        //     , domainParticles{initializer()}
        //     , patchGhostParticles{initializer()}
        //     , levelGhostParticles{initializer()}
        //     , levelGhostParticlesOld{initializer()}
        //     , levelGhostParticlesNew{initializer()}
        //     , pack{name,
        //            &domainParticles,
        //            &patchGhostParticles,
        //            &levelGhostParticles,
        //            &levelGhostParticlesOld,
        //            &levelGhostParticlesNew}
        //     , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
        //     , name_{name}
        // {
        //     validate_ghosts(ghost);
        // }

        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
                      std::string const& name)
            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , domainParticles{construct_particles(ghost[0])}
            , patchGhostParticles{construct_particles(ghost[0])}
            , levelGhostParticles{construct_particles(ghost[0])}
            , levelGhostParticlesOld{construct_particles(ghost[0])}
            , levelGhostParticlesNew{construct_particles(ghost[0])}
            , pack{name,
                   &domainParticles,
                   &patchGhostParticles,
                   &levelGhostParticles,
                   &levelGhostParticlesOld,
                   &levelGhostParticlesNew}
            , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
            , name_{name}
        {
            validate_ghosts(ghost);
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

                using enum core::LayoutMode;
                auto& soa = [&]() -> auto& {
                    if constexpr (any_in(ParticleArray::layout_mode, AoS, AoSTS, AoSMapped, SoATS))
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
                // particles.reserve(n_particles);
                for (std::size_t i = 0; i < n_particles; ++i)
                    particles.emplace_back(soa.copy(i));

                particles.check();
            };

            getParticles("domainParticles", domainParticles);
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
            SamBox const intersectionBox{sourceBox * myGhostBox};

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



        template<typename... Args>
        void copy_from_ghost(Args&&... args);


        void copy_from_cell_overlap(ParticlesData const& pSource,
                                    SAMRAI::pdat::CellOverlap const& pOverlap)
        {
            SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
            SAMRAI::hier::BoxContainer const& boxList = pOverlap.getDestinationBoxContainer();
            for (auto const& overlapBox : boxList)
                copy_(overlapBox, pSource, transformation);
        }

        /**
         * @brief copy with an overlap given by SAMARAI.
         * At runtime we can deal with two kinds of overlaps:
         * - ParticlesDomainOverlap: means this copy is from a context when we're grabbing
         *   leaving domain particles from the neighbor patch, in the patchghost array.
         * - CellOverlap: means domain particles are copied as part of a refinement operation.
         */


        void copy(SAMRAI::hier::PatchData const& source,
                  SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticlesData::copy with overlap");

            // casts throw on failure
            auto& pSource = dynamic_cast<ParticlesData const&>(source);

            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                copy_from_ghost(pSource, *particleOverlap);

            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                copy_from_cell_overlap(pSource, *pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
        }


        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
                   [[maybe_unused]] SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            throw std::runtime_error("Cannot cast");
        }



        bool canEstimateStreamSizeFromBox() const override { return false; }

        std::size_t getOutGoingDataStreamSize(ParticlesDomainOverlap const& pOverlap) const
        {
            auto& transformation        = pOverlap.getTransformation();
            auto const& offset          = as_point<dim>(transformation);
            auto const& noffset         = offset * -1;
            std::size_t numberParticles = 0;
            for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
                numberParticles += patchGhostParticles.nbr_particles_in(
                    shift(phare_box_from<dim>(overlapBox), noffset));
            return sizeof(std::size_t) + numberParticles * ParticleArray::size_of_particle();
        }


        std::size_t getCellOverlapDataStreamSize(SAMRAI::pdat::CellOverlap const& pOverlap) const
        {
            return sizeof(std::size_t)
                   + countNumberParticlesIn_(pOverlap) * ParticleArray::size_of_particle();
        }

        std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                return getOutGoingDataStreamSize(*particleOverlap);

            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                return getCellOverlapDataStreamSize(*pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
        }




        void pack_from_ghost(SAMRAI::tbox::MessageStream&, ParticlesDomainOverlap const&) const;

        void pack_from_cell_overlap(SAMRAI::tbox::MessageStream& stream,
                                    SAMRAI::pdat::CellOverlap const& pOverlap) const
        {
            using PackArray = core::ParticleArray<ParticleArray::options.with_layout(
                core::base_layout_type<ParticleArray>())>;

            if (pOverlap.isOverlapEmpty())
            {
                constexpr std::size_t zero = 0;
                stream << zero;
            }
            else
            {
                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                PackArray outBuffer; // make thread static
                pack_(pOverlap, transformation, outBuffer);
                stream << outBuffer.size();
                stream.growBufferAsNeeded();
                this->pack_(stream, outBuffer);
            }
        }

        /**
         * @brief packStream is the function that takes particles from our particles arrays
         * that lie in the boxes of the given overlap, and pack them to a stream.
         *
         * Streaming particles means that we have to take particles with iCell on a local source
         * index space , communicate them, and load them at destination with iCell on a
         * destination local index space. To do that we need to:
         *
         * 1- translate source iCell to source AMR index space
         * 2- Apply the offset to shift this AMR index on top of the destination cells
         * 3- pack and communicate particles
         * 4- move back iCell from the shifted AMR index space to the local destination index
         * space
         *
         * Note that step 2 could be done upon reception of the pack, we chose to do it before.
         *
         * As for copy(), we can have two kinds of overlaps:
         * - ParticlesDomainOverlap : for grabbing leaving domain particles
         * - CellOverlap : copy as part of refinement operations
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::packStream");

            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
            {
                pack_from_ghost(stream, *particleOverlap);
            }
            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                pack_from_cell_overlap(stream, *pOverlap);
            else
                throw std::runtime_error("Unknown overlap type");
        }


        template<typename ParticleArray_t>
        void pack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t const& outBuffer) const
        {
            if constexpr (any_in(ParticleArray_t::layout_mode, core::LayoutMode::SoA))
                std::apply(
                    [&](auto const&... container) {
                        ((stream.pack(container.data(), outBuffer.size())), ...);
                    },
                    outBuffer.as_tuple());
            else
                stream.pack(outBuffer.data(), outBuffer.size());
        }



        void unpack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                               ParticlesDomainOverlap const& overlap);

        void unpack_cell_overlap(SAMRAI::tbox::MessageStream& stream,
                                 SAMRAI::pdat::CellOverlap const& pOverlap)
        {
            if (pOverlap.isOverlapEmpty())
                return;

            if constexpr (ParticleArray::is_mapped)
                unpackMappedParticles(stream, pOverlap);
            else
                unpackAnyParticles(stream, pOverlap);
        }


        /**
         * @brief unpackStream is the function that unpacks a stream of particles to our
         * domain particle array
         *
         * We get a stream and an overlap. The overlap contains boxes where to put particles and
         * transformation from source to destination AMR indexes.
         *
         * By convention chosen in patckStream, packed particles have their iCell in our AMR
         * index space since we are the destination.
         *
         * like for packStream, we can have two kinds of overlaps:
         * - ParticlesDomainOverlap : for unpacking leaving domain particles
         * - CellOverlap : unpacking as part of refinement operations
         *
         */

        template<typename... Args>
        void unpack_from_ghost(Args&&... args);

        void unpackMappedParticles(SAMRAI::tbox::MessageStream& stream,
                                   SAMRAI::hier::BoxOverlap const& overlap);

        void unpackAnyParticles(SAMRAI::tbox::MessageStream& stream,
                                SAMRAI::hier::BoxOverlap const& overlap);


        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::unpackStream");

            if (auto* particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                unpack_from_ghost(stream, *particleOverlap);

            else if (auto const* pOverlap
                     = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                unpack_cell_overlap(stream, *pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
        }

        template<typename ParticleArray_t>
        void unpack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t& specie) const
        {
            using enum core::LayoutMode;
            if constexpr (any_in(ParticleArray_t::layout_mode, SoA /*, SoATS*/))
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
    };


} // namespace amr
} // namespace PHARE

namespace PHARE::amr
{



template<typename ParticleArray>
template<typename SamraiBox>
void ParticlesData<ParticleArray>::copyPartitionedParticles(
    SamraiBox const& intersectionBox, ParticlesData<ParticleArray> const& sourceData)
{
    auto myDomainBox        = this->getBox();
    auto const per_particle = [&](auto const& particle) {
        if (isInBox(intersectionBox, particle))
            domainParticles.push_back(particle);
    };

    std::array particlesArrays{&sourceData.domainParticles /*, &sourceData.patchGhostParticles*/};

    // for each particles in the source ghost and domain particle arrays
    // we check if it is in the intersectionBox
    // if it is, is it in my domain box ?
    //      - if so, let's add it to my domain particle array
    //      - if not, let's add it to my ghost particle array

    using enum core::LayoutMode;

    for (auto const& sourceParticlesArray : particlesArrays)
    {
        if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS))
        {
            for (auto const& tile : (*sourceParticlesArray)())
                for (auto const& particle : tile())
                    per_particle(particle);
        }
        else
            for (auto const& particle : *sourceParticlesArray)
                per_particle(particle);
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
    auto myDomainBox = this->getBox();
    auto offset      = transformation.getOffset();

    auto const per_particle = [&](auto const& particle) {
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
    };

    std::array const particlesArrays{&sourceData.domainParticles, &sourceData.patchGhostParticles};
    using enum core::LayoutMode;
    for (auto const& sourceParticlesArray : particlesArrays)
    {
        if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS))
        {
            for (auto const& tile : (*sourceParticlesArray)())
                for (auto const& particle : tile())
                    per_particle(particle);
        }
        else
            for (auto const& particle : *sourceParticlesArray)
                per_particle(particle);
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
    if (transformation.getRotation() != SAMRAI::hier::Transformation::NO_ROTATE)
        throw std::runtime_error("Error - rotations not handled in PHARE");

    core::Point const shift{core::for_N<dim, core::for_N_R_mode::make_array>(
        [&](auto i) { return transformation.getOffset()[i]; })};

    auto const offset               = transformation.getOffset();
    auto const per_particle_overlap = [&](auto const& particle, auto const& intersectionBox) {
        auto shiftedParticle = particle.copy();
        for (auto i = 0u; i < dim; ++i)
            shiftedParticle.iCell()[i] += offset[i];
        if (isInBox(intersectionBox, shiftedParticle))
            outBuffer.push_back(shiftedParticle);
    };


    SAMRAI::hier::BoxContainer const& boxContainer = pOverlap.getDestinationBoxContainer();
    auto const& sourceGhostBox                     = getGhostBox();

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

    std::array const particlesArrays{&domainParticles, &patchGhostParticles};
    using enum core::LayoutMode;
    for (auto const& overlapBox : boxContainer)
    {
        SamBox const intersectionBox{transformedSource * overlapBox};

        for (auto const& sourceParticlesArray : particlesArrays)
        {
            if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS))
            {
                for (auto const& tile : (*sourceParticlesArray)())
                    for (auto const& particle : tile())
                        per_particle_overlap(particle, intersectionBox);
            }
            else
                for (auto const& particle : *sourceParticlesArray)
                    per_particle_overlap(particle, intersectionBox);
        }
    }
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

template<typename ParticleArray>
void ParticlesData<ParticleArray>::unpackMappedParticles(SAMRAI::tbox::MessageStream& stream,
                                                         SAMRAI::hier::BoxOverlap const& overlap)
{
    using UnpackArray = core::ParticleArray<ParticleArray::options.with_layout(
        core::base_layout_type<ParticleArray>())>;

    auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

    if (!pOverlap.isOverlapEmpty())
    {
        // unpack particles into a particle array
        std::size_t numberParticles = 0;
        stream >> numberParticles;
        UnpackArray particleArray{numberParticles};
        unpack_(stream, particleArray);

        // ok now our goal is to put the particles we have just unpacked
        // into the particleData and in the proper particleArray : interior or ghost

        SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            // we loop over all boxes in the overlap
            // we have to first take the intersection of each of these boxes
            // with our ghostBox. This is where unpacked particles should go.

            SAMRAI::hier::BoxContainer const& overlapBoxes = pOverlap.getDestinationBoxContainer();

            auto myBox = getBox();

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
void ParticlesData<ParticleArray_t>::unpackAnyParticles(SAMRAI::tbox::MessageStream& stream,
                                                        SAMRAI::hier::BoxOverlap const& overlap)
{
    using UnpackArray = core::ParticleArray<ParticleArray_t::options.with_layout(
        core::base_layout_type<ParticleArray_t>())>;

    auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};
    if (pOverlap.isOverlapEmpty())
        return;

    if (pOverlap.getTransformation().getRotation() != SAMRAI::hier::Transformation::NO_ROTATE)
        return;

    auto const myBox                = getBox();
    auto const per_particle_overlap = [&](auto const& particle, auto const& intersect) {
        if (isInBox(intersect, particle))
        {
            if (isInBox(myBox, particle))
                domainParticles.push_back(particle);
            else
                patchGhostParticles.push_back(particle);
        }
    };

    std::size_t numberParticles = 0;
    stream >> numberParticles;
    UnpackArray particleArray(numberParticles);
    unpack_(stream, particleArray);

    using enum core::LayoutMode;
    SAMRAI::hier::BoxContainer const& overlapBoxes = pOverlap.getDestinationBoxContainer();
    for (auto const& overlapBox : overlapBoxes)
    {
        auto const intersect = getGhostBox() * overlapBox;

        for (auto const& particle : particleArray)
            per_particle_overlap(particle, intersect);
    }
}



template<typename ParticleArray_t>
template<typename... Args>
void ParticlesData<ParticleArray_t>::copy_from_ghost(Args&&... args)
{
    PHARE_LOG_SCOPE(3, "ParticlesData::copy_from_ghost");

    auto&& [pSource, pOverlap] = std::forward_as_tuple(args...);
    auto& src_particles        = pSource.patchGhostParticles;
    auto& dst_particles        = domainParticles;
    auto const& offset         = as_point<dim>(pOverlap.getTransformation());
    auto const& noffset        = offset * -1;

    // auto const offseter = [&](auto const& particle) {
    //     auto shiftedParticle{particle};
    //     for (std::size_t idir = 0; idir < dim; ++idir)
    //         shiftedParticle.iCell()[idir] += offset[idir];
    //     return shiftedParticle;
    // };
    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        core::select_particles(src_particles, dst_particles,
                               shift(phare_box_from<dim>(overlapBox), noffset), offset);
    // src_particles.export_particles(shift(phare_box_from<dim>(overlapBox), noffset),
    //                                dst_particles, offseter);
}




template<typename ParticleArray_t>
void ParticlesData<ParticleArray_t>::pack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                                                     ParticlesDomainOverlap const& pOverlap) const
{
    using PackArray = core::ParticleArray<ParticleArray_t::options.with_layout(
        core::base_layout_type<ParticleArray_t>())>;

    PHARE_LOG_SCOPE(3, "ParticlesData::pack_from_ghost");

    if (pOverlap.isOverlapEmpty())
    {
        constexpr std::size_t zero = 0;
        stream << zero;
        return;
    }

    PackArray outBuffer;
    auto& src_particles = patchGhostParticles;
    auto const& offset  = as_point<dim>(pOverlap.getTransformation());
    auto const& noffset = offset * -1;

    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        core::select_particles(src_particles, outBuffer,
                               shift(phare_box_from<dim>(overlapBox), noffset), offset);

    stream << outBuffer.size();
    stream.growBufferAsNeeded();
    stream.pack(outBuffer.data(), outBuffer.size());
}



// The overlap is not needed here as the pack selects only from the desired overlap
//  and the transform if applicable is performed during packing
template<typename ParticleArray_t>
void ParticlesData<ParticleArray_t>::unpack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                                                       ParticlesDomainOverlap const& /*pOverlap*/)
{
    using PackArray = core::ParticleArray<ParticleArray_t::options.with_layout(
        core::base_layout_type<ParticleArray_t>())>;

    PHARE_LOG_SCOPE(3, "ParticlesData::unpack_from_ghost");

    std::size_t numberParticles = 0;
    stream >> numberParticles;
    PackArray particleArray(numberParticles);
    stream.unpack(particleArray.data(), numberParticles);

    for (auto const& p : particleArray)
        domainParticles.push_back(p);
}



} // namespace PHARE::amr



#endif
