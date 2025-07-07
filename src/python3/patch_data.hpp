#ifndef PHARE_PYTHON_PATCH_DATA_HPP
#define PHARE_PYTHON_PATCH_DATA_HPP

#include <string>
#include <cstring>
#include <utility>

#include "pybind_def.hpp"

namespace PHARE::pydata
{
namespace py = pybind11;

template<typename Data, std::size_t dim>
struct PatchData
{
    static auto constexpr dimension = dim;

    PatchData() = default;

    template<typename... Args>
    PatchData(Args&&... args)
        : data{std::forward<Args...>(args...)}
    {
    }

    Data data;
    std::string patchID;
    py_array_t<std::size_t> origin{dim};
    py_array_t<std::size_t> lower{dim};
    py_array_t<std::size_t> upper{dim};
    std::size_t nGhosts;
};


template<typename PatchData, typename Container>
void setPatchData(PatchData& data, std::string const patchID, Container const origin,
                  Container const lower, Container const upper)
{
    constexpr std::size_t bytes = PatchData::dimension * sizeof(size_t);
    std::memcpy(data.lower.request().ptr, lower.data(), bytes);
    std::memcpy(data.upper.request().ptr, upper.data(), bytes);
    std::memcpy(data.origin.request().ptr, origin.data(), bytes);
    data.patchID = patchID;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string const& patchID)
{
    setPatchData(pdata, patchID, grid.origin().template toArray<std::size_t>(),
                 grid.AMRBox().lower.template toArray<std::size_t>(),
                 grid.AMRBox().upper.template toArray<std::size_t>());
}

template<typename PatchData, typename Field, typename GridLayout>
void setPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                           std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    pdata.data.assign(field.data(), field.data() + field.size());
}


template<typename PatchData, typename Field, typename GridLayout>
void setPyPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                             std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));

    static_assert(PatchData::dimension >= 1 and PatchData::dimension <= 2 && "3d todo");

    if constexpr (PatchData::dimension == 1)
        pdata.data = py::memoryview::from_buffer( //
            field.data(),                         // buffer pointer
            field.shape(),                        // shape (rows, cols)
            {sizeof(double)}                      // strides in bytes
        );

    else if constexpr (PatchData::dimension == 2)
        pdata.data = py::memoryview::from_buffer(               //
            field.data(),                                       // buffer pointer
            field.shape(),                                      // shape (rows, cols)
            {sizeof(double) * field.shape()[1], sizeof(double)} // strides in bytes
        );
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_DATA_H*/
