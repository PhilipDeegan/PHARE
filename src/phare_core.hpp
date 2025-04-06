#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include "core/logger.hpp"
#include "core/vector.hpp"

#include "core/data/grid/grid.hpp"

#include "core/data/electromag/electromag.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/models/physical_state.hpp"
#include "core/models/physical_state.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/algorithm.hpp"

#include "cppdict/include/dict.hpp"


#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <unordered_map>

namespace PHARE::core
{
// = core::LayoutMode::AoSMapped, auto allocator_mode_ = AllocatorMode::CPU

// allows ignoring "interp_order" for types that don't need it
template<std::size_t dim, std::size_t interp_order_ = 0, //
         auto L_ = LayoutMode::AoSMapped, auto A_ = AllocatorMode::CPU>
struct PHARE_Types
{
    static_assert(std::is_same_v<decltype(L_), LayoutMode>);
    static_assert(std::is_same_v<decltype(A_), AllocatorMode>);

    auto static constexpr allocator_mode = A_;
    bool static constexpr c_ordering     = true;
    auto static constexpr dimension      = dim;
    auto static constexpr interp_order   = interp_order_;

    using ParticleSoA_t = SoAParticleArray<dim>;
    using ParticleArray_t
        = ParticleArray<dim, ParticleArrayInternals<dim, L_, StorageMode::VECTOR, A_>>;

    using Array_t          = NdArrayVector<dim, double, c_ordering, allocator_mode>;
    using Grid_t           = Grid<Array_t, HybridQuantity::Scalar>;
    using Field_t          = Field<dim, HybridQuantity::Scalar, double, allocator_mode>;
    using VecField_t       = VecField<Field_t, HybridQuantity>;
    using Electromag_t     = Electromag<VecField_t>;
    using SymTensorField_t = SymTensorField<Field_t, HybridQuantity>;
    using GridLayout_t     = GridLayout<GridLayoutImplYee<dimension, interp_order>>;
    using IonPopulation_t  = IonPopulation<ParticleArray_t, VecField_t, SymTensorField_t>;
    using Ions_t           = Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t      = Electrons<Ions_t>;
    using ParticleInitializerFactory_t = ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
};

struct PHARE_Sim_Types
{
    using SimFunctorParams = cppdict::Dict<int, unsigned int, double, std::size_t>;
    using SimFunctor       = std::function<void(SimFunctorParams const& /*params*/)>;
    using SimulationFunctors // code place id -> function_id -> function
        = std::unordered_map<std::string, std::unordered_map<std::string, SimFunctor>>;
};

} // namespace PHARE::core

#endif // PHARE_CORE_INCLUDE_HPP
