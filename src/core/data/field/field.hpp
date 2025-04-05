#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>

#include "core/def.hpp"

#include "core/data/ndarray/ndarray_view.hpp"

namespace PHARE::core::basic
{
template<std::size_t dim, typename PhysicalQuantity, typename Data_t = double>
class Field : public NdArrayView<dim, Data_t> // NO STRINGS OR STD LIB!
{
public:
    using value_type = Data_t;
    using Super      = NdArrayView<dim, Data_t>;

    Field(PhysicalQuantity qty, value_type* data = nullptr,
          std::array<std::uint32_t, dim> const& dims = ConstArray<std::uint32_t, dim>())
        : Super{data, dims}
        , qty_{qty}
    {
    }

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }

    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }

protected:
    PhysicalQuantity qty_;
};

} // namespace PHARE::core::basic

namespace PHARE::core
{
//! Class Field represents a multidimensional (1,2 or 3D) scalar field
/** Users of Field objects needing to know which physical quantity a specific
 *  Field instance represents can get this info by calling physicalQuantity().
 *  Users may also give a string name to a field object and get a name by calling
 *  name().
 */
template<std::size_t dim, typename PhysicalQuantity, typename Data_t = double,
         auto alloc_mode_ = AllocatorMode::CPU>
class Field : public basic::Field<dim, PhysicalQuantity, Data_t>
{
    static_assert(std::is_same_v<decltype(alloc_mode_), AllocatorMode>);

public:
    using Super                      = basic::Field<dim, PhysicalQuantity, Data_t>;
    auto constexpr static dimension  = dim;
    auto constexpr static alloc_mode = alloc_mode_;
    using value_type                 = Data_t;
    using physical_quantity_type     = PhysicalQuantity;


    Field(std::string const& name, PhysicalQuantity qty, value_type* data = nullptr,
          std::array<std::uint32_t, dim> const& dims = ConstArray<std::uint32_t, dim>())
        : Super{qty, data, dims}
        , name_{name}
    {
    }

    Field(Field const& source)            = default;
    Field(Field&& source)                 = default;
    Field& operator=(Field&& source)      = default;
    Field& operator=(Field const& source) = default;

    auto& operator=(Field* src)
    {
        setBuffer(src);
        return *this;
    }


    NO_DISCARD auto& name() const { return name_; }

    void copyData(Field const& source) { Super::fill_from(source); }



    void setBuffer(std::nullptr_t ptr) _PHARE_ALL_FN_ { setBuffer(static_cast<Field*>(nullptr)); }

    template<typename FieldLike>
    void setBuffer(FieldLike* const field) _PHARE_ALL_FN_
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            assert(field->name() == this->name());
            Super::setShape(field->shape());
        }
        Super::setBuffer(data);
    }

    void setData(Data_t* const data) _PHARE_ALL_FN_ { Super::setBuffer(data); }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }


    template<typename... Args>
    NO_DISCARD auto& operator()(Args&&... args) _PHARE_ALL_FN_
    {
        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            PHARE_DEBUG_DO(                                                                 //
                if (!isUsable()) throw std::runtime_error("Field is not usable: " + name_); //
            )
        }

        return super()(args...);
    }
    template<typename... Args>
    NO_DISCARD auto const& operator()(Args&&... args) const _PHARE_ALL_FN_
    {
        return super()(args...);
    }

    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

private:
    std::string name_{"No Name"};

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }
};



template<std::size_t dim, typename PhysicalQuantity, typename Data_t, auto alloc_mode>
void average(Field<dim, PhysicalQuantity, Data_t, alloc_mode> const& f1,
             Field<dim, PhysicalQuantity, Data_t, alloc_mode> const& f2,
             Field<dim, PhysicalQuantity, Data_t, alloc_mode>& avg)
{
    std::transform(std::begin(f1), std::end(f1), std::begin(f2), std::begin(avg),
                   std::plus<double>());

    std::transform(std::begin(avg), std::end(avg), std::begin(avg),
                   [](double x) { return x * 0.5; });
}


} // namespace PHARE::core


#endif
