#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
// #include "core/data/field/field.hpp"
// #include "core/data/grid/gridlayoutdefs.hpp"
// #include "core/data/particles/particle_array_def.hpp"
#include "core/def/phare_config.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

// #include "core/vector.hpp"
#include "core/utilities/equality.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include <cstdint>

namespace PHARE::core
{

/*
A UsableTensorField is an extension of the TensorField view that owns memory for components and sets
the view pointers. It is useful for tests to easily declare usable (== set views) tensors
*/
template<std::size_t dim, std::size_t rank = 2, auto alloc_mode = AllocatorMode::CPU>
class UsableTensorField : public TensorField<Field_t<dim, alloc_mode>, HybridQuantity, rank>
{
    static_assert(std::is_same_v<decltype(alloc_mode), AllocatorMode>);

    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank>();

    using NdArrayVector_t = NdArrayVector<dim, double, /*c_ordering =*/true, alloc_mode>;
    using Grid_t          = Grid<NdArrayVector_t, HybridQuantity::Scalar>;

public:
    auto static constexpr dimension = dim;
    using Super                     = TensorField<Field_t<dim, alloc_mode>, HybridQuantity, rank>;

protected:
    using tensor_t = typename Super::tensor_t;

    void _set()
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            super()[i].setBuffer(&xyz[i]);
    }


public:
    template<typename GridLayout>
    UsableTensorField(std::string const& name, GridLayout const& layout, tensor_t qty, double v = 0)
        : Super{name, qty}
        , xyz{make_grids(Super::componentNames(), layout, qty)}
    {
        if (v != 0)
            for (auto& grid : xyz)
                grid.fill(v);
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
            tensorfield[i].setBuffer(&xyz[i]);
    }


    Super& super() { return *this; }
    Super& super() const { return *this; }

    Super& view() { return *this; }
    Super const& view() const { *this; }

    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }

    auto& grids() const { return xyz; }

    template<auto AM>
    bool isclose(UsableTensorField<dim, rank, AM> const& that, double diff = 1e-15) const
    {
        return core::for_N_all<N_elements>(
            [&](auto i) { return (*this)[i].isclose(that[i], diff); });
    }

    template<auto AM>
    bool operator==(UsableTensorField<dim, rank, AM> const& that) const
    {
        return core::for_N_all<N_elements>([&](auto i) { return (*this)[i] == that[i]; });
    }


protected:
    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = HybridQuantity::componentsQuantities(qty);
        return for_N<N_elements, for_N_R_mode::make_array>(
            [&](auto i) { return Grid_t{compNames[i], qts[i], layout.allocSize(qts[i])}; });
    }

    std::array<Grid_t, N_elements> xyz;
};


template<bool binary_eq = false, typename F0, typename F1, typename PhysicalQuantity,
         std::size_t rank>
EqualityReport compare_tensor_fields(TensorField<F0, PhysicalQuantity, rank> const& ref,
                                     TensorField<F1, PhysicalQuantity, rank> const& cmp,
                                     [[maybe_unused]] double const diff)
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank>();

    auto const same_sizes = [&]() {
        return core::for_N_all<N_elements>([&](auto i) { return ref[i].size() == cmp[i].size(); });
    }();

    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    auto const float_eq = [&](auto const a, auto const b) {
        if constexpr (binary_eq)
            return a == b;
        else
            return float_equals(a, b, diff);
    };

    std::stringstream log;

    for (std::size_t ci = 0; ci < N_elements; ++ci)
    {
        auto const& ref_c                       = ref[ci];
        auto const& cmp_c                       = cmp[ci];
        auto const [eqvals, eqnot0, ref0, cmp0] = FloatFieldComparator_t{diff}(ref_c, cmp_c);
        if (eqvals != ref_c.size())
        {
            auto const bad = ref_c.size() - eqvals;
            log << "Tensorfield value mismatch: \n";
            log << " component: " << ci << " - ";
            log << "ok(" << eqvals << ") - ";
            log << "ok!=0(" << eqnot0 << ") - ";
            log << "bad(" << bad << ") - ";
            log << "ref0(" << ref0 << ") - ";
            log << "cmp0(" << cmp0 << ")\n";
            return EqualityReport{false, log.str()};
        }
    }

    return EqualityReport{true};
}


template<bool binary_equal = false, typename Field_t, typename PhysicalQuantity, std::size_t rank>
EqualityReport compare_tensor_fields_strict(TensorField<Field_t, PhysicalQuantity, rank> const& ref,
                                            TensorField<Field_t, PhysicalQuantity, rank> const& cmp,
                                            double const diff)
{
    if (ref.componentNames() != cmp.componentNames())
        return EqualityReport{false, "Tensorfield component mismatch"};

    // else check rest
    return compare_tensor_fields(ref, cmp, diff);
}


} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
