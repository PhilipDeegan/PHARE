#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"

#include <array>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace PHARE::core
{
template<std::size_t dim, bool c_ordering = true>
struct NdArrayViewer
{
    using Id  = std::uint16_t;
    using Idx = Id const;

    template<typename NCells, template<typename, std::size_t> typename Indexes, typename Index>
    static inline std::uint32_t idx(NCells const& nCells,
                                    Indexes<Index, dim> const& indexes) _PHARE_ALL_FN_
    {
        if constexpr (dim == 1)
            return idx(nCells, indexes[0]);

        else if constexpr (dim == 2)
            return idx(nCells, indexes[0], indexes[1]);

        else if constexpr (dim == 3)
            return idx(nCells, indexes[0], indexes[1], indexes[2]);
    }

    static inline std::uint32_t idx(auto const /*nCells*/, Idx i) _PHARE_ALL_FN_ { return i; }


    static inline std::uint32_t idx(auto const nCells, Idx i, Idx j) _PHARE_ALL_FN_
    {
        if constexpr (c_ordering)
            return j + i * nCells[1];
        else
            return i + j * nCells[0];
    }
    static inline std::uint32_t idx(auto const nCells, Idx i, Idx j, Idx k) _PHARE_ALL_FN_
    {
        if constexpr (c_ordering)
            return k + j * nCells[2] + i * nCells[1] * nCells[2];
        else
            return i + j * nCells[0] + k * nCells[1] * nCells[0];
    }



    template<template<typename, std::size_t> typename Indexes, typename Index>
    NO_DISCARD static inline auto& at(auto* data, auto const& nCells,
                                      Indexes<Index, dim> const& indexes) _PHARE_ALL_FN_

    {
        auto const& i = idx(nCells, indexes);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }

    static inline auto& at(auto* data, auto const nCells, auto const... indexes) _PHARE_ALL_FN_
    {
        auto const& i = idx(nCells, indexes...);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }
};



template<std::size_t dim, typename DataType = double, bool c_ordering = true>
class NdArrayView
{
    using viewer = NdArrayViewer<dim, c_ordering>;

public:
    std::size_t static const dimension = dim;
    using type                         = DataType;
    using value_type                   = DataType;
    using pointer_type                 = DataType*;


    NdArrayView(pointer_type ptr, std::array<std::uint32_t, dim> const nCells) _PHARE_ALL_FN_
        : ptr_{ptr},
          size_{core::product(nCells)},
          nCells_{nCells}
    {
    }


    NdArrayView(NdArrayView const&)            = default;
    NdArrayView& operator=(NdArrayView const&) = default;


    template<typename Index>
    NO_DISCARD inline auto& operator[](std::array<Index, dim> const& indexes) _PHARE_ALL_FN_
    {
        return viewer::at(ptr_, nCells_, indexes);
    }
    template<typename Index>
    NO_DISCARD inline auto& operator[](std::array<Index, dim> const& indexes) const _PHARE_ALL_FN_
    {
        return viewer::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD inline auto const&
    operator()(std::array<Index, dim> const& indexes) const _PHARE_ALL_FN_
    {
        return viewer::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD inline auto& operator()(std::array<Index, dim> const& indexes) _PHARE_ALL_FN_
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes));
    }


    NO_DISCARD inline auto const& operator()(auto const... indexes) const _PHARE_ALL_FN_
    {
        return viewer::at(ptr_, nCells_, indexes...);
    }

    inline auto& operator()(auto const... indexes) _PHARE_ALL_FN_
    {
        return viewer::at(ptr_, nCells_, indexes...);
    }


    NO_DISCARD auto& data() const _PHARE_ALL_FN_ { return ptr_; }
    NO_DISCARD auto& data() _PHARE_ALL_FN_ { return ptr_; }

    NO_DISCARD auto& size() const _PHARE_ALL_FN_ { return size_; }
    NO_DISCARD auto& shape() const _PHARE_ALL_FN_ { return nCells_; }

    NO_DISCARD auto begin() const _PHARE_ALL_FN_ { return ptr_; }
    NO_DISCARD auto begin() _PHARE_ALL_FN_ { return ptr_; }

    NO_DISCARD auto end() const _PHARE_ALL_FN_ { return ptr_ + size_; }
    NO_DISCARD auto end() _PHARE_ALL_FN_ { return ptr_ + size_; }


    void zero() _PHARE_ALL_FN_ { fill(1e-36); }
    auto zeros() const
    {
        return sum_from(*this, [](auto const e) { return e == 0 ? 1 : 0; });
    }
    // void check() const
    // {
    //     PHARE_DEBUG_DO({
    //         for (auto const& e : *this)
    //         {
    //             if (std::isnan(e))
    //                 throw std::runtime_error("NAN");
    //             if (std::isinf(e))
    //                 throw std::runtime_error("INF");
    //         }
    //     })
    // }

    // void notZero() const
    // {
    //     for (std::size_t i = 0; i < size(); ++i)
    //         if (std::abs(data()[i]) < 1e-15)
    //             throw std::runtime_error("ZERO");
    // }

    auto& fill(DataType const& v) _PHARE_ALL_FN_
    {
        std::fill(begin(), end(), v);
        return *this;
    }

    bool isclose(NdArrayView const& that, double diff = 1e-13) const
    {
        if (this->size() != that.size())
            return false;
        for (std::size_t i = 0; i < this->size(); ++i)
            if (!float_equals(this->data()[i], that.data()[i], diff))
                return false;
        return true;
    }


    bool operator==(NdArrayView const& that) const
    {
        if (this->size() != that.size())
            return false;
        for (std::size_t i = 0; i < this->size(); ++i)
            if (this->data()[i] != that.data()[i])
                return false;
        return true;
    }

    template<typename View>
    void reset(View& view) _PHARE_ALL_FN_
    {
        this->ptr_    = view.data();
        this->size_   = view.size();
        this->nCells_ = view.nCells_;
    }
    template<typename Vec>
    void reset(Vec& vec, std::array<std::uint32_t, dim> const& nCells) _PHARE_ALL_FN_
    {
        this->ptr_    = vec.data();
        this->size_   = vec.size();
        this->nCells_ = nCells;
    }

    void fill_from(NdArrayView const& that)
    {
        if (for_N_any<dim>([&](auto i) { return shape()[i] != that.shape()[i]; }))
            throw std::runtime_error("ArrayView::fill_from: Incompatible input shape");
        std::copy(that.data(), that.data() + size(), data());
    }

    void setBuffer(pointer_type ptr) _PHARE_ALL_FN_ { ptr_ = ptr; }
    void setShape(std::array<std::uint32_t, dim> const nCells) _PHARE_ALL_FN_
    {
        nCells_ = nCells;
        size_   = core::product(nCells);
    }

    auto& reshape(auto const& shape)
    {
        setShape(shape);
        return *this;
    }


    auto size_address() _PHARE_ALL_FN_ { return &size_; }

private:
    pointer_type ptr_ = nullptr;
    std::size_t size_;
    std::array<std::uint32_t, dim> nCells_;
};


template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType* data, std::array<std::uint32_t, dim> const shape) _PHARE_ALL_FN_
{
    return NdArrayView<dim, DataType, c_ordering>{data, shape};
}

template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType const* const data,
                     std::array<std::uint32_t, dim> const shape) _PHARE_ALL_FN_
{
    return NdArrayView<dim, DataType const, c_ordering>{data, shape};
}

template<typename DataType, std::size_t dim>
auto make_array_view(std::vector<DataType>& vec, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType>{vec.data(), shape};
}


template<typename DataType, std::size_t dim>
auto make_array_view(std::vector<DataType> const& vec, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType const>{vec.data(), shape};
}




} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP
