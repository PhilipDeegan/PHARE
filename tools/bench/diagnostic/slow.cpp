#include <memory>
#define PHARE_FORCE_LOG_LINE 1


#include "core/logger.hpp"
#include "core/utilities/mpi_utils.hpp"

#include <tuple>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include <stdexcept>

#include "highfive/highfive.hpp"


std::string const FILE_NAME("slooooow.h5");
std::string const DATASET_NAME("dset");

/* float byte size table
     1MB = 262120
    10MB = 2621200
    50MB = 13106000
    75MB = 19659000
   100MB = 26212000
   200MB = 52424000
   500MB = 131060000
*/

std::array static inline const SIZES = {
    262120ull, 2621200ull, 13106000ull, 19659000ull, 26212000ull, 52424000ull, 131060000ull,
};

namespace PHARE::core
{

std::uint64_t static now()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

struct Timer
{
    ~Timer()
    {
        if (mpi::rank() == 0)
            PHARE_LOG_LINE_SS(id << " : time : " << now() - start);
    }

    auto operator()() const { return now() - start; }

    std::string id;
    std::size_t start = now();
};

auto rank_bits(auto const size)
{
    std::size_t per_rank = size / mpi::size();
    std::size_t offset   = per_rank * mpi::rank();
    return std::make_tuple(per_rank, offset);
}

struct CSV
{
    std::ofstream static setup(auto const& csv)
    {
        bool const newFile = !std::ofstream(csv, std::ios::out | std::ios::in);

        std::ofstream of(csv, std::ios::out | std::ios::app);
        if (newFile)
            of << "time" << std::endl;
        of.flush();
        return of;
    }


    void operator<<(auto const line) { of << line << std::endl; }

    std::string const csv = "csv." + std::to_string(mpi::rank()) + ".csv";
    std::ofstream of{setup(csv)};
};



// clang-format off
#define TIME(s) Timer PHARE_STR_CAT(_, __LINE__){s}
// clang-format on

void run(auto const size)
{
    using namespace HighFive;

    if (size % mpi::size() != 0)
        throw std::runtime_error("size is not divisible by ranks");

    CSV csv;
    TIME("tools/bench/diagnostic/slow.cpp       ");


    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
    std::unique_ptr<File> fileptr
        = std::make_unique<File>(FILE_NAME, File::ReadWrite | File::Create | File::Truncate, fapl);
    auto& file = *fileptr;

    // precreate
    file.createDataSet(DATASET_NAME, DataSpace({size}, {size}), create_datatype<float>());

    Timer timer{"tools/bench/diagnostic/slow.cpp::write"};
    {
        auto dataset                    = file.getDataSet(DATASET_NAME);
        auto const&& [per_rank, offset] = rank_bits(size);
        std::vector<float> data(per_rank, offset);
        for (std::size_t i = 0; i < per_rank; ++i)
            data[i] += i;
        dataset.select({offset}, {per_rank}).write(data);
        csv << (std::to_string(size) + "," + std::to_string(timer()));
    }
    fileptr.reset();

    if (mpi::rank() == 0)
        CSV{"total.csv"} << (std::to_string(size) + "," + std::to_string(timer()));
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    for (auto const size : SIZES)
        PHARE::core::run(size);
    MPI_Finalize();

    return 0;
}
