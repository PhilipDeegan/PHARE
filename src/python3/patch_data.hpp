#ifndef PHARE_PYTHON_PATCH_DATA_HPP
#define PHARE_PYTHON_PATCH_DATA_HPP


#include "pybind_def.hpp"
#include "core/utilities/point/point.hpp"


#include <string>
#include <cstring>
#include <utility>


namespace PHARE::pydata
{
namespace py = pybind11;

template<typename Data, std::size_t dim>
struct __attribute__((visibility("hidden"))) PatchData
{
    static auto constexpr dimension = dim;

    PatchData() = default;

    template<typename... Args>
    PatchData(Args&&... args)
        requires std::is_constructible_v<Data, Args&&...>
        : data{std::forward<Args>(args)...}
    {
    }

    Data data;
    std::string patchID;
    py_array_t<double> origin{dim};
    py_array_t<int> lower{dim};
    py_array_t<int> upper{dim};
    std::size_t nGhosts;
};


template<typename PatchData>
void setPatchData(PatchData& data, std::string const patchID, auto const origin, auto const lower,
                  auto const upper)
{
    std::memcpy(data.lower.request().ptr, lower.data(), PatchData::dimension * sizeof(int));
    std::memcpy(data.upper.request().ptr, upper.data(), PatchData::dimension * sizeof(int));
    std::memcpy(data.origin.request().ptr, origin.data(), PatchData::dimension * sizeof(double));
    data.patchID = patchID;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string const& patchID)
{
    setPatchData(pdata, patchID, *grid.origin(), *grid.AMRBox().lower, *grid.AMRBox().upper);
}


template<typename PatchData, typename Field, typename GridLayout>
void setPyPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                             std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));

    static_assert(PatchData::dimension >= 1 and PatchData::dimension <= 3);

    pdata.data = py::memoryview::from_buffer(field.data(), field.shape(), strides(field.shape()));
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_DATA_H*/
