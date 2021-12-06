
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_range.h"
#include "test_ranges.h"
#include "test_range_replacer.h"



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
