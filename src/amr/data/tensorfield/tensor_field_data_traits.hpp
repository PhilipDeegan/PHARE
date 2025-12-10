#ifndef PHARE_SRC_AMR_TENSOR_FIELD_TENSOR_FIELD_DATA_TRAITS_HPP
#define PHARE_SRC_AMR_TENSOR_FIELD_TENSOR_FIELD_DATA_TRAITS_HPP

#include <array>
#include <concepts>

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Patch.h"

namespace PHARE
{
namespace amr
{
    /**
     * @brief Concept ensuring a type satisfies the PHARE FieldData interface.
     * This allows templates like RefineFieldPatchStrategy to accept either
     * FieldData or TensorFieldData.
     */
    template<typename T>
    concept IsTensorFieldData
        = std::derived_from<T, SAMRAI::hier::PatchData>
          && requires(T a, T const ca, SAMRAI::hier::Patch const& patch) {
                 // Type aliases
                 typename T::gridlayout_type;
                 typename T::grid_type;
                 typename T::tensor_field_type;

                 // Static constexpr variables
                 requires std::same_as<decltype(T::dimension), std::size_t const>;
                 requires std::same_as<decltype(T::interp_order), std::size_t const>;
                 requires std::same_as<decltype(T::N), std::size_t const>;

                 // Public member variables
                 requires std::same_as<decltype(a.gridLayout), typename T::gridlayout_type>;
                 requires std::same_as<decltype(a.grids), std::array<typename T::grid_type, T::N>>;

                 // API requirements
                 { a.getPointer() } -> std::same_as<std::array<typename T::grid_type, T::N>*>;
                 { T::getLayout(patch, 0) } -> std::same_as<typename T::gridlayout_type const&>;
                 {
                     T::getFields(patch, 0)
                 } -> std::same_as<std::array<typename T::grid_type, T::N>&>;
                 { T::getTensorField(patch, 0) } -> std::same_as<typename T::tensor_field_type>;
             };

    template<typename T>
    concept IsVecFieldData = IsTensorFieldData<T> && (T::N == 3);

    /**
     * @brief Compile-time utility to select the correct Geometry type based on field data rank.
     * This meta-function extracts the geometric representation of a field data. It distinguishes
     * between scalar field data and tensor field data to resolve nested type aliases correctly.
     * @tparam ScalarOrTensorFieldDataT The data structure representing the field data.
     * @tparam is_scalar Boolean flag; true if the field is a scalar, false if it is a tensor.
     */
    template<typename ScalarOrTensorFieldDataT, bool is_scalar>
    struct FieldGeometrySelector;
    /**
     * @brief Specialization for scalar field data
     */
    template<typename ScalarOrTensorFieldDataT>
    struct FieldGeometrySelector<ScalarOrTensorFieldDataT, true>
    {
        using type = ScalarOrTensorFieldDataT::Geometry;
    };
    /**
     * @brief Specialization for tensor field data
     */
    template<typename ScalarOrTensorFieldDataT>
    struct FieldGeometrySelector<ScalarOrTensorFieldDataT, false>
    {
        using type = ScalarOrTensorFieldDataT::Geometry::FieldGeometry_t;
    };


} // namespace amr
} // namespace PHARE

#endif // PHARE_SRC_AMR_TENSOR_FIELD_TENSOR_FIELD_DATA_TRAITS_HPP
