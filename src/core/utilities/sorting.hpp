#ifndef PHARE_SORTING_HPP
#define PHARE_SORTING_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


#include <iostream>
#include <cstddef>
#include <vector>

namespace PHARE::core
{

template<typename Array, std::size_t dim>
class CountingSort
{
public:
    CountingSort(std::size_t size = 0)
        : tmp_(size)
    {
    }

    using value_type = typename Array::value_type;
    using SortBox    = core::Box<int, dim>;

    void setup(std::size_t nbr_elements, SortBox const& box)
    {
        tmp_.resize(nbr_elements);
        hist_.resize(box.size());
        zero(hist_);
    }

    template<typename CellGetter>
    void sort(Array& toSort, CellGetter getCell)
    {
        // ABORT_IF_NOT(toSort.size() <= tmp_.size());

        // compute histogram
        for (std::size_t ip = 0; ip < tmp_.size(); ++ip)
        {
            auto const& item = toSort[ip];
            auto const& cell = getCell(item);
            // if (cell > hist_.size())
            //     PHARE_LOG_LINE_STR(Point{item.iCell} << " " << toSort.box());
            // ABORT_IF(cell > hist_.size());
            hist_[cell]++;
        }

        std::uint32_t sum = 0;
        for (std::size_t i = 0; i < hist_.size(); ++i)
        {
            // all particles in cell i will be in [sum, sum+hist_[i])
            int const tmp = hist_[i];
            hist_[i]      = sum;
            sum += tmp;
        }

        for (std::size_t ip = 0; ip < toSort.size(); ++ip)
        {
            auto const& cell    = getCell(toSort[ip]);
            tmp_[hist_[cell]++] = toSort[ip];
        }

        // now put items back in toSort
        // in the right order
        for (std::size_t ip = 0; ip < toSort.size(); ++ip)
        {
            toSort[ip] = tmp_[ip];
        }
    }

    // used below to convert to ndarray
    auto& histogram() const { return hist_; }

private:
    std::vector<value_type> tmp_;
    std::vector<std::uint32_t> hist_;
};



template<typename Array_t>
struct CountingSortCellConverter
{
    auto constexpr static dim = Array_t::dimension;
    using NdArray_t           = NdArrayVector<dim, std::uint32_t>;
    using box_t               = core::Box<int, dim>;

    template<typename CellGetter>
    void operator()(box_t const& box, CellGetter getCell)
    {
        auto const& histogram = sorter.histogram();
        // ABORT_IF_NOT(histogram.size() == box.size());

        std::size_t i = 0, prev = 0;
        for (auto const& cell : box)
        {
            auto n_items                = histogram[i] - prev;
            particles.ppc(cell)         = n_items;
            particles.ppc_offsets(cell) = histogram[i] - n_items;
            prev                        = histogram[i];
            ++i;
        }
    }

    CountingSort<Array_t, dim>& sorter;
    Array_t& particles;
};


} // namespace PHARE::core



#endif
