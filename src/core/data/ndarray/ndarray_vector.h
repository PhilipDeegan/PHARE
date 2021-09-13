#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H


#include <iostream>

#include <stdexcept>
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <numeric>

#include "core/def.h"

#ifdef HAVE_UMPIRE
#include "umpire/ResourceManager.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"
#include "kul/log.hpp"
#endif

namespace PHARE::core
{
template<std::size_t dim, typename DataType = double>
struct NdArrayViewer
{
    template<typename NCells, typename... Indexes>
    static DataType const& at(DataType const* data, NCells const& nCells,
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




template<std::size_t dim, typename DataType = double, typename Pointer = DataType const*>
class NdArrayView : NdArrayViewer<dim, DataType>
{
public:
    static constexpr bool is_contiguous = 1;
    static const std::size_t dimension  = dim;
    using type                          = DataType;

    NdArrayView(Pointer ptr, std::array<std::uint32_t, dim> const nCells) _PHARE_ALL_FN_
        : ptr_{ptr},
          nCells_{nCells}
    {
    }

    NdArrayView(std::vector<DataType> const& v, std::array<std::uint32_t, dim> const& nbCell)
        : NdArrayView{v.data(), nbCell}
    {
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes... indexes) const _PHARE_ALL_FN_
    {
        return NdArrayViewer<dim, DataType>::at(ptr_, nCells_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes... indexes) _PHARE_ALL_FN_
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes...));
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

    auto data() const { return ptr_; }
    auto data() { return ptr_; }

    std::size_t size() const
    {
        return std::accumulate(nCells_.begin(), nCells_.end(), 1, std::multiplies<std::size_t>());
    }
    auto shape() const { return nCells_; }

private:
    Pointer ptr_ = nullptr;
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
            assert(rm.isAllocator("samrai::data_allocator"));
            return Allocator{rm.getAllocator("samrai::data_allocator")};
        }
#endif
    }

} // namespace

template<std::size_t dim, typename DataType = double,
         typename Allocator_ = typename std::vector<DataType>::allocator_type>
class NdArrayVector
{
    static std::uint32_t accumulate(std::array<std::uint32_t, dim> const& ncells)
    {
        return std::accumulate(ncells.begin(), ncells.end(), 1, std::multiplies<std::uint32_t>());
    }

    void _print()
    {
        KLOG(INF) << data() << " " << size();
    }

public:
    static constexpr bool is_contiguous = 1;
    static constexpr bool is_host_mem   = is_host_mem_type<DataType, Allocator_>();
    static const std::size_t dimension  = dim;
    using Allocator                     = Allocator_;
    using type                          = DataType;

    NdArrayVector() = delete;

    template<typename... Nodes>
    explicit NdArrayVector(Nodes... nodes)
        : size_{(... * nodes)}
        , nCells_{nodes...}
        , data_(get_allocator<DataType, Allocator>())
    {
        data_.resize(size_);
        static_assert(sizeof...(Nodes) == dim);
        _print();
    }

    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells)
        : size_{accumulate(ncells)}
        , nCells_{ncells}
        , data_(get_allocator<DataType, Allocator>())
    {
        data_.resize(size_);
        _print();
    }

    NdArrayVector(NdArrayVector const& source) = default;
    NdArrayVector(NdArrayVector&& source)      = default;

    auto size() const { return size_; }

    auto data() { return data_.data(); }
    auto data() const { return data_.data(); }

    auto begin() const { return std::begin(data_); }
    auto begin() { return std::begin(data_); }

    auto end() const { return std::end(data_); }
    auto end() { return std::end(data_); }

    void zero()
    {
        if(size() == 0) return;
        if constexpr (is_host_mem)
            std::fill(data_.begin(), data_.end(), 0);
#if defined(HAVE_RAJA)
        else
        {
            KLOG(INF) << size();
            std::vector<DataType> zeroes(size(), 0);
            assert(zeroes.size() == size());
            // TODO replace with "fill"
            RAJA::resources::Cuda{}.memcpy(
                /*device pointer*/ data(),
                /*host pointer*/ zeroes.data(),
                /*size in bytes*/ sizeof(DataType) * size());
        }
#endif
    }


    NdArrayVector& operator=(NdArrayVector const& source)
    {
        if (nCells_ != source.nCells_)
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");

        if constexpr (is_host_mem)
            this->data_ = source.data_;
        else
            assert(false);

        return *this;
    }

    NdArrayVector& operator=(NdArrayVector&& source)
    {
        if (nCells_ != source.nCells_)
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");

        if constexpr (is_host_mem)
            this->data_ = std::move(source.data_);
        else
            assert(false);

        return *this;
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(data(), nCells_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes...));
    }

    template<typename Index>
    DataType const& operator()(std::array<Index, dim> const& indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(data(), nCells_, indexes);
    }

    template<typename Index>
    DataType& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes));
    }


    auto shape() const { return nCells_; }

    template<typename Mask>
    auto operator[](Mask&& mask)
    {
        return MaskedView{*this, std::forward<Mask>(mask)};
    }

private:
    std::size_t size_;
    std::array<std::uint32_t, dim> nCells_;

    std::vector<DataType, Allocator> data_;
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
