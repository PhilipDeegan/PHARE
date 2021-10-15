#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H

#include <array>
#include <tuple>
#include <vector>
#include <cstdint>
#include <numeric>
#include <stdexcept>

#include "core/def.h"
#include "core/utilities/types.h"

#ifdef HAVE_UMPIRE
#include "umpire/ResourceManager.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"
#include "kul/log.hpp"
#include "SAMRAI/hier/ForAll.h"
#endif

namespace PHARE::core
{
template<std::size_t dim, typename DataType = double>
struct NdArrayViewer
{
    template<typename NCells, typename... Indexes>
    static auto& at(DataType const* data, NCells const& nCells,
                    Indexes const&... indexes) _PHARE_ALL_FN_
    {
        auto params = std::forward_as_tuple(indexes...);
        static_assert(sizeof...(Indexes) == dim);
        // static_assert((... && std::is_unsigned_v<decltype(indexes)>)); TODO : manage later if
        // this test should be included

        if constexpr (dim == 1)
        {
            auto i = std::get<0>(params);

            return data[i];
        }

        if constexpr (dim == 2)
        {
            auto i = std::get<0>(params);
            auto j = std::get<1>(params);

            return data[j + i * nCells[1]];
        }

        if constexpr (dim == 3)
        {
            auto i = std::get<0>(params);
            auto j = std::get<1>(params);
            auto k = std::get<2>(params);

            return data[k + j * nCells[2] + i * nCells[1] * nCells[2]];
        }
    }

    template<typename NCells, typename Index>
    static DataType const& at(DataType const* data, NCells const& nCells,
                              std::array<Index, dim> const& indexes) _PHARE_ALL_FN_

    {
        if constexpr (dim == 1)
            return data[indexes[0]];

        else if constexpr (dim == 2)
            return data[indexes[1] + indexes[0] * nCells[1]];

        else if constexpr (dim == 3)
            return data[indexes[2] + indexes[1] * nCells[2] + indexes[0] * nCells[1] * nCells[2]];
    }
};



template<typename Array, typename Mask>
class MaskedView
{
public:
    static auto constexpr dimension = Array::dimension;
    using DataType                  = typename Array::type;
    using data_type                 = typename Array::type;

    MaskedView(Array& array, Mask const& mask)
        : array_{array}
        , shape_{array.shape()}
        , mask_{mask}
    {
    }

    MaskedView(Array& array, Mask&& mask)
        : array_{array}
        , shape_{array.shape()}
        , mask_{std::move(mask)}
    {
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dimension, DataType>::at(array_.data(), shape_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<MaskedView const&>(*this)(indexes...));
    }

    auto operator=(data_type value) { mask_.fill(array_, value); }

    auto xstart() const { return mask_.min(); }

    auto xend() const { return shape_[0] - 1 - mask_.max(); }


    auto ystart() const { return mask_.min(); }

    auto yend() const { return shape_[1] - 1 - mask_.max(); }


private:
    Array& array_;
    std::array<std::uint32_t, dimension> shape_;
    Mask const& mask_;
};




template<std::size_t dim, typename DataType = double, typename Ptr_ = DataType* const,
         bool is_host_mem_ = true>
class NdArrayView
{
public:
    static constexpr auto is_contiguous = true;
    static constexpr auto is_host_mem   = is_host_mem_;
    static const std::size_t dimension  = dim;
    using type                          = DataType;
    using pointer_type                  = Ptr_;
    using This                          = NdArrayView<dim, DataType, Ptr_>;
    using view_t                        = This;

    NdArrayView(pointer_type ptr, std::array<std::uint32_t, dim> const nCells) _PHARE_ALL_FN_
        : ptr_{ptr},
          size_{core::product(nCells)},
          nCells_{nCells}
    {
    }


    template<typename Index>
    DataType const& operator()(std::array<Index, dim> const& indexes) const _PHARE_ALL_FN_
    {
        return NdArrayViewer<dim, DataType>::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    DataType& operator()(std::array<Index, dim> const& indexes) _PHARE_ALL_FN_
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes));
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes const... indexes) const _PHARE_ALL_FN_
    {
        return NdArrayViewer<dim, DataType>::at(ptr_, nCells_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes const... indexes) _PHARE_ALL_FN_
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes...));
    }


    auto& data() const _PHARE_ALL_FN_ { return ptr_; }
    auto& data() _PHARE_ALL_FN_ { return ptr_; }

    auto& size() const _PHARE_ALL_FN_ { return size_; }
    auto& shape() const _PHARE_ALL_FN_ { return nCells_; }

    auto begin() const _PHARE_ALL_FN_ { return ptr_; }
    auto begin() _PHARE_ALL_FN_ { return ptr_; }

    auto end() const _PHARE_ALL_FN_ { return ptr_ + size_; }
    auto end() _PHARE_ALL_FN_ { return ptr_ + size_; }

private:
    pointer_type ptr_ = nullptr;
    std::size_t size_;
    std::array<std::uint32_t, dim> nCells_;
};

namespace
{
    template<typename DataType, typename Allocator>
    constexpr bool is_host_mem_type()
    {
        if constexpr (std::is_same_v<Allocator, typename std::vector<DataType>::allocator_type>)
            return true;
        return false;
    }


    template<typename DataType, typename Allocator>
    auto get_allocator()
    {
        auto constexpr is_host_mem = is_host_mem_type<DataType, Allocator>();

        if constexpr (is_host_mem)
            return Allocator{};
#if defined(HAVE_UMPIRE)
        else
        {
            auto& rm = umpire::ResourceManager::getInstance();
            assert(rm.isAllocator("PHARE::data_allocator"));
            return Allocator{rm.getAllocator("PHARE::data_allocator")};
        }
#endif
    }

    template<typename vector_impl, typename Allocator>
    auto initialize(std::size_t size, Allocator&& allocator)
    {
        vector_impl vec(std::forward<Allocator>(allocator));
        vec.resize(size);
        return vec;
    }

} // namespace

template<std::size_t dim, typename DataType = double,
         typename Allocator_ = typename std::vector<DataType>::allocator_type>
class NdArrayVector
    : public StackVar<std::vector<DataType>>,
      public NdArrayView<dim, DataType, DataType* const, is_host_mem_type<DataType, Allocator_>()>
{
public:
    static constexpr bool is_contiguous = 1;
    static const std::size_t dimension  = dim;

    using Allocator   = Allocator_;
    using vector_impl = std::vector<DataType, Allocator>;
    using Vector      = StackVar<vector_impl>;
    using Super       = NdArrayView<dim, DataType>;
    using type        = DataType;
    using view_t      = Super;

    using Super::data;
    using Super::is_host_mem;
    using Super::shape;
    using Super::size;
    using Vector::var;

    template<typename... Nodes>
    explicit NdArrayVector(Nodes... nodes)
        : Vector{initialize<vector_impl>((... * nodes), get_allocator<DataType, Allocator>())}
        , Super{Vector::var.data(), {nodes...}}
    {
        static_assert(sizeof...(Nodes) == dim);
    }

    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells)
        : Vector{initialize<vector_impl>(core::product(ncells),
                                         get_allocator<DataType, Allocator>())}
        , Super{Vector::var.data(), ncells}
    {
    }


    NdArrayVector(NdArrayVector const& source) = default;
    NdArrayVector(NdArrayVector&& source)      = default;



    NdArrayVector& operator=(NdArrayVector const& source)
    {
        if (shape() != source.shape())
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");

        this->var = source.var;
        return *this;
    }

    NdArrayVector& operator=(NdArrayVector&& source)
    {
        if (shape() != source.shape())
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");

        this->var = std::move(source.var);
        return *this;
    }


    Super const& view() const { return *this; }
    Super& view() { return *this; }


    template<typename Mask>
    auto operator[](Mask&& mask)
    {
        return MaskedView{*this, std::forward<Mask>(mask)};
    }

    void zero()
    {
        if (size() == 0)
            return;
        if constexpr (is_host_mem)
            std::fill(this->var.begin(), this->var.end(), 0);
#if defined(HAVE_RAJA)
        else
            RAJA::resources::Cuda{}.memset(data(), 0, sizeof(DataType) * size());
#endif
    }
};


class NdArrayMask
{
public:
    NdArrayMask(std::size_t min, std::size_t max)
        : min_{min}
        , max_{max}
    {
    }

    NdArrayMask(std::size_t width)
        : min_{width}
        , max_{width}
    {
    }

    template<typename Array>
    void fill(Array& array, typename Array::type val) const
    {
        if constexpr (Array::dimension == 1)
            fill1D(array, val);

        else if constexpr (Array::dimension == 2)
            fill2D(array, val);

        else if constexpr (Array::dimension == 3)
            fill3D(array, val);
    }

    template<typename Array>
    void fill1D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        for (std::size_t i = min_; i <= max_; ++i)
            array(i) = val;

        for (std::size_t i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            array(i) = val;
    }

    template<typename Array>
    void fill2D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        // left border
        for (std::size_t i = min_; i <= max_; ++i)
            for (std::size_t j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;

        // right border
        for (std::size_t i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            for (std::size_t j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;


        for (std::size_t i = min_; i <= shape[0] - 1 - min_; ++i)
        {
            // bottom border
            for (std::size_t j = min_; j <= max_; ++j)
                array(i, j) = val;

            // top border
            for (std::size_t j = shape[1] - 1 - max_; j <= shape[1] - 1 - min_; ++j)
                array(i, j) = val;
        }
    }

    template<typename Array>
    void fill3D(Array& array, typename Array::type val) const
    {
        throw std::runtime_error("3d not implemented");
    }

    template<typename Array>
    auto nCells(Array const& array)
    {
        auto shape = array.shape();

        std::size_t cells = 0;

        if constexpr (Array::dimension == 1)
            for (std::size_t i = min_; i <= max_; ++i)
                cells += 2;

        if constexpr (Array::dimension == 2)
            for (std::size_t i = min_; i <= max_; ++i)
                cells += (shape[0] - (i * 2) - 2) * 2 + (shape[1] - (i * 2) - 2) * 2 + 4;

        if constexpr (Array::dimension == 3)
            throw std::runtime_error("Not implemented dimension");

        return cells;
    }


    auto min() const { return min_; };
    auto max() const { return max_; };

private:
    std::size_t min_, max_;
};




template<typename Array, typename Mask>
void operator>>(MaskedView<Array, Mask>&& inner, MaskedView<Array, Mask>&& outer)
{
    using MaskedView_t = MaskedView<Array, Mask>;

    if constexpr (MaskedView_t::dimension == 1)
    {
        assert(inner.xstart() > outer.xstart());
        assert(inner.xend() < outer.xend());
        outer(outer.xstart()) = inner(inner.xstart());
        outer(outer.xend())   = inner(inner.xend());
    }


    if constexpr (MaskedView_t::dimension == 2)
    {
        assert(inner.xstart() > outer.xstart() and inner.xend() < outer.xend()
               and inner.ystart() > outer.ystart() and inner.yend() < outer.yend());

        for (auto ix = inner.xstart(); ix <= inner.xend(); ++ix)
        {
            outer(ix, outer.ystart()) = inner(ix, inner.ystart()); // bottom
            outer(ix, outer.yend())   = inner(ix, inner.yend());   // top
        }

        for (auto iy = inner.ystart(); iy <= inner.yend(); ++iy)
        {
            outer(outer.xstart(), iy) = inner(inner.xstart(), iy); // left
            outer(outer.xend(), iy)   = inner(inner.xend(), iy);   // right
        }

        // bottom left
        for (auto ix = outer.xstart(); ix < inner.xstart(); ++ix)
            outer(ix, outer.ystart()) = inner(inner.xstart(), inner.ystart());

        for (std::size_t iy = outer.ystart(); iy < inner.ystart(); ++iy)
            outer(outer.xstart(), iy) = inner(inner.xstart(), inner.ystart());


        // top left
        for (auto ix = outer.xstart(); ix < inner.xstart(); ++ix)
            outer(ix, outer.yend()) = inner(inner.xstart(), inner.yend());

        for (auto iy = outer.yend(); iy > inner.yend(); --iy)
            outer(outer.xstart(), iy) = inner(inner.xstart(), inner.yend());

        // top right
        for (auto ix = outer.xend(); ix > inner.xend(); --ix)
            outer(ix, outer.yend()) = inner(inner.xend(), inner.yend());

        for (auto iy = outer.yend(); iy > inner.yend(); --iy)
            outer(outer.xend(), iy) = inner(inner.xend(), inner.yend());


        // bottom right
        for (auto ix = outer.xend(); ix > inner.xend(); --ix)
            outer(ix, outer.ystart()) = inner(inner.xend(), inner.ystart());

        for (auto iy = outer.ystart(); iy < inner.ystart(); ++iy)
            outer(outer.xend(), iy) = inner(inner.xend(), inner.ystart());
    }

    if constexpr (MaskedView_t::dimension == 3)
    {
        throw std::runtime_error("3d not implemented");
    }
}

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
