#ifndef PHARE_CORE_DATA_FIELD_FIELD_TRAITS
#define PHARE_CORE_DATA_FIELD_FIELD_TRAITS

#include <concepts>
#include <cstddef>
#include <string>
#include <type_traits>

namespace PHARE::core
{

template<typename T>
concept IsField = requires(T field) {
    // 1. Static Metadata
    { T::dimension } -> std::convertible_to<std::size_t>;
    typename T::value_type;
    typename T::physical_quantity_type;

    // 2. Identification Interface
    { field.name() } -> std::convertible_to<std::string const&>;
    { field.physicalQuantity() } -> std::same_as<typename T::physical_quantity_type const&>;

    // 3. State/Memory Management
    { field.isUsable() } -> std::same_as<bool>;
    { field.data() } -> std::same_as<typename T::value_type*>; // Inherited from NdArrayView

    // 4. Access Pattern (checking for 1D, 2D, or 3D access capability)
    // We use a requires expression to ensure it accepts 'dimension' number of arguments
    requires((T::dimension == 1 && requires(T f) {
                 { f(std::declval<std::size_t>()) } -> std::same_as<typename T::value_type&>;
             }) || (T::dimension == 2 && requires(T f) {
                 {
                     f(std::declval<std::size_t>(), std::declval<std::size_t>())
                 } -> std::same_as<typename T::value_type&>;
             }) || (T::dimension == 3 && requires(T f) {
                 {
                     f(std::declval<std::size_t>(), std::declval<std::size_t>(),
                       std::declval<std::size_t>())
                 } -> std::same_as<typename T::value_type&>;
             }));
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_FIELD_FIELD_TRAITS
