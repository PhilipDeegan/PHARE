#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP


#include <array>
#include <tuple>
#include <vector>
#include <cstdint>
#include <numeric>
#include <iostream>
#include <stdexcept>


#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/types.hpp"


#include "core/data/ndarray/ndarray_view.hpp"
#include "core/data/ndarray/ndarray_mask.hpp"


namespace PHARE::core
{

template<typename Type, auto allocator_mode_>
struct NdArrayConsts
{
    auto static constexpr allocator_mode = allocator_mode_;

    using vector_t  = typename Vector<Type, allocator_mode_>::vector_t;
    using stack_var = StackVar<vector_t>;
};


template<std::size_t dim, typename Type = double, bool c_ordering = true,
         auto allocator_mode_ = AllocatorMode::CPU,
         typename Consts      = NdArrayConsts<Type, allocator_mode_>>
class NdArrayVector : public Consts::stack_var, public NdArrayView<dim, Type>
{
public:
    auto static constexpr allocator_mode = allocator_mode_;
    auto static constexpr dimension      = dim;

    using Storage = typename Consts::stack_var;
    using View    = NdArrayView<dim, Type>;
    using Storage::var;
    using value_type = Type;
    using View::data;
    using View::shape;
    using View::size;

    using vec_helper = PHARE::Vector<Type, allocator_mode>;

    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells)
        : Storage{vec_helper::make(core::product(ncells))}
        , View{Storage::var.data(), ncells}
    {
    }


    NdArrayVector(NdArrayVector const& that)
        : Storage{vec_helper::from(that.var)}
        , View{Storage::var.data(), that.shape()}
    {
    }


    template<std::size_t d, typename T, bool c, auto a, typename C>
    NdArrayVector(NdArrayVector<d, T, c, a, C>&& that)
        : Storage{vec_helper::from(std::move(that.var))}
        , View{Storage::var.data(), that.shape()}
    {
    }


    template<typename... Nodes,
             typename = typename std::enable_if_t<(std::is_integral_v<std::decay_t<Nodes>> && ...)>>
    explicit NdArrayVector(Nodes&&... nodes)
        : NdArrayVector{std::array<std::uint32_t, dim>{nodes...}}
    {
        static_assert(sizeof...(Nodes) == dim);
    }

    auto& operator=(std::array<std::uint32_t, dim> const& ncells)
    {
        Storage::var = std::move(vec_helper::make(core::product(ncells)));
        View::reset(this->var, ncells);
        return *this;
    }


    auto& operator=(NdArrayVector const& that)
    {
        if constexpr (allocator_mode == AllocatorMode::CPU)
            this->var = that.var;
        else
            vec_helper::copy(this->var, that.var);

        View::reset(this->var, that.shape());

        return *this;
    }

    auto& operator=(NdArrayVector&& that)
    {
        if constexpr (allocator_mode == AllocatorMode::CPU)
            this->var = std::move(that.var);
        else
            vec_helper::copy(this->var, that.var);

        View::reset(this->var, that.shape());

        return *this;
    }



    NO_DISCARD auto operator[](NdArrayMask&& mask)
    {
        return MaskedView{*this, std::forward<NdArrayMask>(mask)};
    }
    NO_DISCARD auto operator[](NdArrayMask const& mask) { return MaskedView{*this, mask}; }



    void zero()
    {
        if (size() == 0)
            return;
        vec_helper::fill(this->var, 0);
    }


    auto& vector() { return Storage::var; }
    auto& vector() const { return Storage::var; }

    void reset() { View::reset(this->var, shape()); }

    View const& view() const { return *this; }
    View& view() { return *this; }

    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }
};




template<typename F, std::size_t dim, typename Type, bool c_ordering, auto alloc_mode>
auto& update_from(F f, NdArrayVector<dim, Type, c_ordering, alloc_mode> const& in)
{
    for (std::size_t i = 0; i < in.size(); ++i)
        in.data()[i] = f(i);
    return in;
}



template<auto alloc_mode0, typename F, std::size_t dim, typename Type, bool c_ordering,
         auto alloc_mode1>
auto generate_from(F f, NdArrayVector<dim, Type, c_ordering, alloc_mode1> const& in)
{
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    NdArrayVector<dim, value_type, c_ordering, alloc_mode0> ret{in.shape()};
    for (std::size_t i = 0; i < in.size(); ++i)
        ret.data()[i] = f(i);
    return ret;
}


} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP
