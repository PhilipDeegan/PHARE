
#ifndef PHARE_CORE_INCLUDE_H
#define PHARE_CORE_INCLUDE_H

#include "core/data/electromag/electromag.h"
#include "core/data/electrons/electrons.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "core/data/ndarray/ndarray_vector.h"
#include "core/data/particles/particle_array.h"
#include "core/data/vecfield/vecfield.h"
#include "core/models/physical_state.h"
#include "core/models/physical_state.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/algorithm.h"
#include "core/logger.h"

#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <unordered_map>

#include "cppdict/include/dict.hpp"

#if defined(WITH_RAJA) and defined(WITH_UMPIRE)
#include "core/data/particles/llnl/particle_array.h"
#elif defined(PHARE_WITH_GPU)

#endif

namespace PHARE::core
{
template<typename Particle>
auto constexpr particle_array_selector()
{
#if defined(WITH_RAJA) and defined(WITH_UMPIRE)
    return static_cast<PHARE::core::llnl::ParticleArray<Particle>*>(nullptr);
#elif defined(WITH_RAJA) or defined(WITH_UMPIRE)
#error // invalid, both RAJA and UMPIRE are required together.
#else
    return static_cast<PHARE::core::ParticleArray<Particle>*>(nullptr);
#endif
}


template<std::size_t dimension_, std::size_t interp_order_>
struct PHARE_Types
{
    static auto constexpr dimension    = dimension_;
    static auto constexpr interp_order = interp_order_;

    using Array_t      = PHARE::core::NdArrayVector<dimension>;
    using VecField_t   = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t      = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t = PHARE::core::Electromag<VecField_t>;
    using YeeLayout_t  = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;

    using Particle_t      = PHARE::core::Particle<dimension>;
    using ParticleArray_t = std::decay_t<decltype(*particle_array_selector<Particle_t>())>;

    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = PHARE::core::Electrons<Ions_t>;

    using ParticleInitializerFactory
        = PHARE::core::ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
};

struct PHARE_Sim_Types
{
    using SimFunctorParams = cppdict::Dict<int, unsigned int, double, std::size_t>;
    using SimFunctor       = std::function<void(SimFunctorParams const& /*params*/)>;
    using SimulationFunctors // code place id -> function_id -> function
        = std::unordered_map<std::string, std::unordered_map<std::string, SimFunctor>>;
};

} // namespace PHARE::core

#endif // PHARE_CORE_INCLUDE_H
