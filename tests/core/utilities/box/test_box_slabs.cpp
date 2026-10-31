

#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_slabs.hpp"
#include "core/data/field/field_box_spans.hpp"

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <cstdint>
#include <algorithm>
#include <type_traits>


namespace PHARE::core
{



TEST(BoxSlabsTest, test_reverse_iterator)
{
    static constexpr std::size_t dim = 3;
    static constexpr std::size_t IDX = dim - 1;
    static constexpr PHARE::SimOpts opts{dim, 1};
    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;

    GridLayout_t layout{3};

    std::vector<Point<int, dim>> points;

    auto const slabs = make_box_slabs(layout.AMRBox());

    for (auto slabit = slabs.end(); slabit-- > slabs.begin();)
    {
        auto const& slab = *slabit;
        for (auto spans = slab.end(); spans-- > slab.begin();)
        {
            auto const&& [start, size] = *spans;
            auto point                 = start;
            point[IDX] += size;
            for (; point[IDX]-- > start[IDX];)
                points.emplace_back(point);
        }
    }

    std::vector<Point<int, dim>> const expected{
        {2, 2, 2}, {2, 2, 1}, {2, 2, 0}, //
        {2, 1, 2}, {2, 1, 1}, {2, 1, 0}, //
        {2, 0, 2}, {2, 0, 1}, {2, 0, 0}, //
        {1, 2, 2}, {1, 2, 1}, {1, 2, 0}, //
        {1, 1, 2}, {1, 1, 1}, {1, 1, 0}, //
        {1, 0, 2}, {1, 0, 1}, {1, 0, 0}, //
        {0, 2, 2}, {0, 2, 1}, {0, 2, 0}, //
        {0, 1, 2}, {0, 1, 1}, {0, 1, 0}, //
        {0, 0, 2}, {0, 0, 1}, {0, 0, 0},
    };

    EXPECT_EQ(expected, points);
}


TEST(BoxSlabsTest, test_reverse_iterator_unsigned_box)
{
    static constexpr std::size_t dim = 3;
    static constexpr std::size_t IDX = dim - 1;

    Box<std::uint32_t, dim> const box{{0, 0, 0}, {2, 2, 2}};

    std::vector<Point<std::uint32_t, dim>> points;

    auto const slabs = make_box_slabs(box);

    for (auto slabit = slabs.end(); slabit-- > slabs.begin();)
    {
        auto const& slab = *slabit;
        for (auto spans = slab.end(); spans-- > slab.begin();)
        {
            auto const&& [start, size] = *spans;
            auto point                 = start;
            point[IDX] += size;
            for (; point[IDX]-- > start[IDX];)
                points.emplace_back(point);
        }
    }

    std::vector<Point<std::uint32_t, dim>> const expected{
        {2u, 2u, 2u}, {2u, 2u, 1u}, {2u, 2u, 0u}, //
        {2u, 1u, 2u}, {2u, 1u, 1u}, {2u, 1u, 0u}, //
        {2u, 0u, 2u}, {2u, 0u, 1u}, {2u, 0u, 0u}, //
        {1u, 2u, 2u}, {1u, 2u, 1u}, {1u, 2u, 0u}, //
        {1u, 1u, 2u}, {1u, 1u, 1u}, {1u, 1u, 0u}, //
        {1u, 0u, 2u}, {1u, 0u, 1u}, {1u, 0u, 0u}, //
        {0u, 2u, 2u}, {0u, 2u, 1u}, {0u, 2u, 0u}, //
        {0u, 1u, 2u}, {0u, 1u, 1u}, {0u, 1u, 0u}, //
        {0u, 0u, 2u}, {0u, 0u, 1u}, {0u, 0u, 0u},
    };

    EXPECT_EQ(expected, points);
}


TEST(BoxSlabsTest, test_range_loop)
{
    std::size_t static constexpr dim = 3;
    Box<std::uint32_t, dim> box{{0, 0, 0}, {9, 9, 9}};
    std::size_t elements = 0;

    for (auto const& slab : make_box_slabs(box))
        for (auto const& [start, size] : slab)
            elements += size;

    EXPECT_EQ(elements, 10 * 10 * 10);
}



TEST(BoxSlabsTest, test_iter_loop_dim1)
{
    std::size_t static constexpr dim = 1;
    Box<std::uint32_t, dim> box{{0}, {9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_slabs(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;

            elements += size;
        }
    }

    EXPECT_EQ(elements, 10);
}

TEST(BoxSlabsTest, test_iter_loop_dim2)
{
    std::size_t static constexpr dim = 2;
    Box<std::uint32_t, dim> box{{0, 0}, {9, 9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_slabs(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;
            elements += size;
        }
    }

    EXPECT_EQ(elements, 10 * 10);
}

TEST(BoxSlabsTest, test_iter_loop_dim3)
{
    std::size_t static constexpr dim = 3;
    Box<std::uint32_t, dim> box{{0, 0, 0}, {9, 9, 9}};
    std::size_t elements = 0;

    auto const& slabs = make_box_slabs(box);

    for (auto slabit = slabs.begin(); slabit != slabs.end(); ++slabit)
    {
        auto const& slab = *slabit;

        for (auto rowit = slab.begin(); rowit != slab.end(); ++rowit)
        {
            auto const& [start, size] = *rowit;
            elements += size;
        }
    }

    EXPECT_EQ(elements, 10 * 10 * 10);
}


// ill-formed boxes (lower > upper, e.g. an empty SAMRAI intersection) iterate nothing
TEST(BoxSlabsTest, test_ill_formed_box_is_empty)
{
    auto const count = [](auto const& box) {
        std::size_t elements = 0;
        for (auto const& slab : make_box_slabs(box))
            for (auto const& [start, size] : slab)
                elements += size;
        return elements;
    };

    EXPECT_EQ(count(Box<int, 1>{{5}, {3}}), 0);
    EXPECT_EQ(count(Box<int, 2>{{0, 5}, {9, 3}}), 0);
    EXPECT_EQ(count(Box<int, 2>{{5, 0}, {3, 9}}), 0);
    EXPECT_EQ(count(Box<int, 3>{{5, 0, 0}, {3, 9, 9}}), 0);
    EXPECT_EQ(count(Box<int, 3>{{0, 5, 0}, {9, 3, 9}}), 0);
    EXPECT_EQ(count(Box<int, 3>{{0, 0, 5}, {9, 9, 3}}), 0);
    EXPECT_EQ(count(Box<std::uint32_t, 3>{{0, 0, 5}, {9, 9, 3}}), 0);
}


// for_each_point_zip pairs each point of box0, in C order, with the point at the same offset
//  in box1, for boxes of different point types (e.g. a GridLayout's AMR and local boxes)
TEST(BoxSlabsTest, test_for_each_point_zip)
{
    auto const check = [](auto const& box0, auto const& box1) {
        using Point0 = std::decay_t<decltype(box0.lower)>;
        using Point1 = std::decay_t<decltype(box1.lower)>;
        std::vector<Point0> points0;
        std::vector<Point1> points1;

        for_each_point_zip(box0, box1, [&](auto const& p0, auto const& p1) {
            points0.emplace_back(p0);
            points1.emplace_back(p1);
        });

        ASSERT_EQ(points0.size(), box0.size());
        ASSERT_EQ(points1.size(), box1.size());

        // strictly increasing in C order, so every point of box0 is visited once
        auto const not_increasing
            = [](auto const& a, auto const& b) { return !(a.toArray() < b.toArray()); };
        EXPECT_EQ(std::adjacent_find(points0.begin(), points0.end(), not_increasing),
                  points0.end());

        for (std::size_t i = 0; i < points0.size(); ++i)
            for (std::size_t d = 0; d < box0.dimension; ++d)
                EXPECT_EQ(static_cast<long>(points0[i][d]) - static_cast<long>(box0.lower[d]),
                          static_cast<long>(points1[i][d]) - static_cast<long>(box1.lower[d]));
    };

    check(Box<int, 1>{{5}, {9}}, Box<std::uint32_t, 1>{{2}, {6}});
    check(Box<int, 2>{{5, -2}, {9, 1}}, Box<std::uint32_t, 2>{{2, 0}, {6, 3}});
    check(Box<int, 3>{{5, -2, 0}, {7, 1, 2}}, Box<std::uint32_t, 3>{{2, 0, 1}, {4, 3, 3}});
}


// for_each_point_reverse visits exactly the points of for_each_point, last to first, including
//  unsigned boxes starting at 0 (no underflow) and ill-formed boxes (nothing visited)
TEST(BoxSlabsTest, test_for_each_point_reverse)
{
    auto const check = [](auto const& box) {
        using Point_t = std::decay_t<decltype(box.lower)>;
        std::vector<Point_t> forward, reverse;

        for_each_point(box, [&](auto const& p) { forward.emplace_back(p); });
        for_each_point_reverse(box, [&](auto const& p) { reverse.emplace_back(p); });

        std::reverse(forward.begin(), forward.end());
        EXPECT_EQ(forward, reverse);
        return reverse.size();
    };

    EXPECT_EQ(check(Box<int, 1>{{-2}, {9}}), 12);
    EXPECT_EQ(check(Box<int, 2>{{5, -2}, {9, 1}}), 5 * 4);
    EXPECT_EQ(check(Box<int, 3>{{5, -2, 0}, {7, 1, 2}}), 3 * 4 * 3);
    EXPECT_EQ(check(Box<std::uint32_t, 1>{{0}, {4}}), 5);
    EXPECT_EQ(check(Box<std::uint32_t, 2>{{0, 0}, {3, 2}}), 4 * 3);
    EXPECT_EQ(check(Box<std::uint32_t, 3>{{0, 0, 0}, {2, 3, 1}}), 3 * 4 * 2);
    EXPECT_EQ(check(Box<int, 3>{{0, 5, 0}, {9, 3, 9}}), 0);
    EXPECT_EQ(check(Box<std::uint32_t, 2>{{0, 5}, {9, 3}}), 0);
}


// copies of a spans iterator must not share state: advancing/dereferencing one copy leaves
//  the others untouched (see PR #1153 review)
TEST(BoxSlabsTest, test_copies_do_not_alias)
{
    std::size_t static constexpr dim = 3;
    Box<std::uint32_t, dim> const box{{0, 0, 0}, {2, 2, 2}};

    auto const slabs = make_box_slabs(box);
    auto const slab  = *slabs.begin();
    auto const base  = slab.begin();

    auto a = base; // copies of a const iterator
    auto b = base;
    auto c = b; // copy of a non-const iterator

    auto const [a_point, a_size] = *a;
    Point<std::uint32_t, dim> const a_start{a_point};

    ++b;
    auto const [b_point, b_size] = *b;

    EXPECT_NE(&a_point, &b_point);
    EXPECT_EQ(a_point, a_start);
    EXPECT_EQ(a_start, (Point<std::uint32_t, dim>{0u, 0u, 0u}));
    EXPECT_EQ(b_point, (Point<std::uint32_t, dim>{0u, 1u, 0u}));
    EXPECT_EQ(a_size, b_size);

    auto const [c_point, c_size] = *c;
    EXPECT_EQ(c_point, a_start);
    EXPECT_EQ(a_point, a_start);
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
