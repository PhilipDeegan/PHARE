#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_HPP
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_HPP

#include <array>
#include <utility>
#include <algorithm>
#include <unordered_map>

#include "vecfield_component.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

namespace PHARE
{
namespace core
{
    template<typename Field_t, typename PhysicalQuantity>
    using VecField = TensorField<Field_t, PhysicalQuantity, /*rank=*/1>;


    struct VecFieldNames
    {
        std::string vecName;
        std::string xName;
        std::string yName;
        std::string zName;

        VecFieldNames() = default;

        template<typename VecFieldT>
        VecFieldNames(VecFieldT const& v)
            : vecName{v.name()}
            , xName{v.getComponentName(core::Component::X)}
            , yName{v.getComponentName(core::Component::Y)}
            , zName{v.getComponentName(core::Component::Z)}

        {
        }
    };


} // namespace core
} // namespace PHARE

#endif
