
#include "python3/interop.hpp"
#include "python3/interpreter.hpp"

#include "gtest/gtest.h"

TEST(PythonInterop, array_works)
{
    std::size_t constexpr static dim = 2;
    PHARE::pydata::py_array_t<std::size_t> upper{dim};
}


TEST(PythonInterop, Interpreter_can_import_job_file)
{
    PHARE::py3::Interpreter::INSTANCE().import();
}


int main(int argc, char** argv)
{
    PHARE::py3::Interpreter::INSTANCE(); // setup python

    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
