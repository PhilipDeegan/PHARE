#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include "core/data/particles/particle_array_def.hpp"
#include "core/logger.hpp"
#include "core/vector.hpp"

#include "core/data/grid/grid.hpp"

#include "core/data/grid/grid_tiles.hpp"
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
template<typename GridLayout_t, auto L_, auto A_>
struct UsingResolver
{
    auto static constexpr dimension    = GridLayout_t::dimension;
    auto static constexpr interp_order = GridLayout_t::interp_order;

    using Array_t = NdArrayVector<dimension, double, /*c_order*/ true, A_>;
    using Grid_t  = Grid<Array_t, HybridQuantity::Scalar>;
    using Field_t = Field<dimension, HybridQuantity::Scalar, double, A_>;
};

template<typename GridLayout_t, auto A_>
struct UsingResolver<GridLayout_t, LayoutMode::AoSTS, A_>
{
    auto static constexpr dimension    = GridLayout_t::dimension;
    auto static constexpr interp_order = GridLayout_t::interp_order;

    template<typename T, auto am>
    using nd_array_t = NdArrayVector<dimension, T, /*c_order*/ true, am>;
    // template<typename T, auto am>
    // using nd_array_vt = NdArrayView<dimension, T, /*c_order*/ true>;

    using Grid_t  = GridTileSet<GridLayout_t, nd_array_t<double, A_>, HybridQuantity::Scalar>;
    using Field_t = FieldTileSet<GridLayout_t, nd_array_t<double, A_>, HybridQuantity::Scalar>;
};

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

    using GridLayout_t = GridLayout<GridLayoutImplYee<dimension, interp_order>>;

    using Resolver_t = UsingResolver<GridLayout_t, L_, A_>;

    using Array_t = NdArrayVector<dim, double, c_ordering, allocator_mode>;

    using Grid_t  = Resolver_t::Grid_t;  // layout dependent
    using Field_t = Resolver_t::Field_t; // layout dependent

    using VecField_t       = VecField<Field_t, HybridQuantity>;
    using Electromag_t     = Electromag<VecField_t>;
    using SymTensorField_t = SymTensorField<Field_t, HybridQuantity>;
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
