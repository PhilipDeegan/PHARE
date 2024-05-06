#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "phare_core.hpp"
#include "core/data/field/field.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"


namespace PHARE::core
{

template<std::size_t dim, std::size_t rank = 2>
class UsableTensorField : public TensorField<Field_t<dim>, HybridQuantity, rank>
{
public:
    using Super = TensorField<Field_t<dim>, HybridQuantity, rank>;

protected:
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank>();
    using tensor_t                   = typename Super::tensor_t;

public:
    auto static constexpr dimension = dim;

    template<typename GridLayout>
    UsableTensorField(std::string const& name, GridLayout const& layout, tensor_t qty, double v = 0)
        : Super{name, qty}
        , xyz{make_fields(Super::componentNames(), layout, qty)}
    {
        if (v != 0)
            for (auto& grid : xyz)
                std::fill_n(grid->data(), grid->size(), v);
        _set();
    }


    UsableTensorField(UsableTensorField&& that)
        : Super{std::forward<Super>(that)}
        , xyz{std::move(that.xyz)}
    {
        _set();
    }


    void set_on(Super& tensorfield)
    {
        // used for setting on normal model tensorfields
        for (std::size_t i = 0; i < N_elements; ++i)
            tensorfield.setBuffer(Super::componentNames()[i], xyz[i].get());
    }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }
    auto& fields() const { return Super::components(); }

protected:
    template<typename Names, typename GridLayout, std::size_t... Is>
    auto static _make_fields(Names const& compNames, GridLayout const& layout, tensor_t qty,
                             std::integer_sequence<std::size_t, Is...>)
    {
        auto qts = HybridQuantity::componentsQuantities(qty);
        auto fn  = [&](auto i) {
            return std::make_unique<Field_t<dim>>(compNames[i], qts[i], layout.allocSize(qts[i]));
        };
        return std::array<std::unique_ptr<Field_t<dim>>, N_elements>{fn(Is)...};
    }

    template<typename ComponentNames, typename GridLayout>
    auto static make_fields(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        return _make_fields(compNames, layout, qty,
                            std::make_integer_sequence<std::size_t, N_elements>{});
    }

    void _set()
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            Super::setBuffer(Super::componentNames()[i], xyz[i].get());
    }

    std::array<std::unique_ptr<Field_t<dim>>, N_elements> xyz;
};


} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
