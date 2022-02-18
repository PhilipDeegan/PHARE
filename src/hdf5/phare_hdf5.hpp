#ifndef PHARE_HDF5_PHARE_HDF5_HPP
#define PHARE_HDF5_PHARE_HDF5_HPP

#include <tuple>
#include <memory>
#include <vector>

#if PHARE_HAS_HIGHFIVE
#include "highfive/H5Version.hpp"

#define _PHARE_WITH_HIGHFIVE(...) __VA_ARGS__
#else
#define _PHARE_WITH_HIGHFIVE(...)
#endif


#include "core/utilities/types.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


namespace PHARE::hdf5
{
using namespace PHARE::core;

template<typename GridLayout_, typename Model_>
class HybridPatchView
{
    using This = HybridPatchView<GridLayout_, Model_>;

public:
    using GridLayout = GridLayout_;
    using Model      = Model_;

    static constexpr auto dim = GridLayout::dimension;

    using ParticleArray_t = ParticleArray<dim> const&;
    using Field_t         = NdArrayView<dim>;
    using VecField_t      = tuple_fixed_type<NdArrayView<dim>, 3>;

    using IonPopulation_t = std::tuple<VecField_t, Field_t, tuple_fixed_type<ParticleArray_t, 3>>;


    std::unique_ptr<This> static make_unique(Model const& model)
    {
        return std::make_unique<aggregate_adapter<This>>( //
            vecfield_view(model.state.electromag.E),      //
            vecfield_view(model.state.electromag.B),      //
            vecfield_view(model.state.J),                 //
            field_view(model.state.ions.density()),       //
            core::generate(
                [](auto& pop) -> IonPopulation_t {
                    return {vecfield_view(pop.flux()), field_view(pop.density()),
                            std::make_tuple(pop.domainParticles(), pop.patchGhostParticles(),
                                            pop.levelGhostParticles())};
                },
                model.state.ions));
    }


    template<typename Hierarchy>
    auto static levels(Hierarchy const& hierarchy, Model const& model)
    {
        using PatchView = hdf5::HybridPatchView<GridLayout, Model>;

        std::vector<std::vector<std::unique_ptr<PatchView>>> levels(1);

        PHARE::amr::visitHierarchy<GridLayout>(
            hierarchy, *model.resourcesManager,
            [&](auto& layout, auto, auto levelNumber) {
                if (levelNumber >= levels.size())
                    levels.emplace_back();
                auto& level = levels[levelNumber];

                level.emplace_back(PatchView::make_unique(model));
            },
            0, hierarchy.getNumberOfLevels(), model);

        return levels;
    }


    VecField_t const E, B, J;
    Field_t const rho;
    std::vector<IonPopulation_t> const ions;

protected:
    template<typename Field>
    auto static field_view(Field const& field)
    {
        return make_array_view(field);
    }

    template<typename VecField>
    auto static vecfield_view(VecField const& vecField)
    {
        auto const& [x, y, z] = vecField();

        return std::make_tuple(field_view(x), field_view(y), field_view(z));
    }
};


} // namespace PHARE::hdf5


#endif // PHARE_HDF5_PHARE_HDF5_HPP
