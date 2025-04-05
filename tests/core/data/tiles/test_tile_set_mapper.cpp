
#include "core/data/tiles/tile_set_mapper.hpp"

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "test_tile.hpp"
#include <gtest/gtest.h>


template<typename TileSet_>
struct TileSetInfo
{
    auto static constexpr dimension = TileSet_::dimension;
    using Tile_t                    = TileSet_::value_type;
    using Tiles                     = std::vector<Tile_t>;
    using TileSet_t                 = TileSet_;

    Box<int, TileSet_t::dimension> box;
    std::array<std::uint32_t, TileSet_t::dimension> shape;
    std::array<std::size_t, TileSet_t::dimension> tile_size;
};

template<typename TileSet_>
class TileMappingTest : public ::testing::Test
{
public:
    auto static constexpr dimension = TileSet_::dimension;
    using Tile_t                    = TileSet_::value_type;
    using Tiles                     = std::vector<Tile_t>;
    using TileSet_t                 = TileSet_;
    using info_t                    = TileSetInfo<TileSet_t>;

    TileMappingTest() {}

    auto static count_cells(Tiles const& tiles)
    {
        return sum_from(tiles, [](auto const& e) { return e.size(); });
    }

    auto static all_in(Tiles const& tiles, Box<int, dimension> const& box)
    {
        for (auto const& tile : tiles)
            if (box * tile != tile)
                return false;
        return true;
    };


    auto static print(Tiles const& tiles)
    {
        for (auto const& tile : tiles)
        {
            PHARE_LOG_LINE_SS(tile);
        }
    }
};


using DimTiles
    = testing::Types</*TileSet<Box<int, 1>>, TileSet<Box<int, 2>>,*/ TileSet<Box<int, 3>>>;

TYPED_TEST_SUITE(TileMappingTest, DimTiles);

TYPED_TEST(TileMappingTest, simpleSymmetric)
{
    using TileSet_t     = TestFixture::TileSet_t;
    using TileSetInfo_t = TestFixture::info_t;
    auto const& [box, shape, tile_size]
        = TileSetInfo_t{{{0, 0, 0}, {9, 9, 9}}, {2, 2, 2}, {5, 5, 5}};

    for_N<2>([&](auto i) {
        auto const tiles = Tiler<TileSet_t, i>{{box, shape, tile_size}}.f();
        EXPECT_TRUE(this->all_in(tiles, box));
        EXPECT_EQ(this->count_cells(tiles), std::pow(10, this->dimension));
        EXPECT_EQ(tiles.size(), std::pow(2, this->dimension));
    });
}

TYPED_TEST(TileMappingTest, asymmetricImpl0)
{
    using TileSet_t     = TestFixture::TileSet_t;
    using TileSetInfo_t = TestFixture::info_t;
    auto const& [box, shape, tile_size]
        = TileSetInfo_t{{{0, 0, 0}, {18, 16, 14}}, {4, 4, 3}, {5, 5, 5}};

    auto const tiles = Tiler<TileSet_t, 0>{{box, shape, tile_size}}.f();

    EXPECT_TRUE(this->all_in(tiles, box));
    EXPECT_FALSE(any_overlaps_in(tiles, [](auto const& tile) { return tile; }));
    EXPECT_EQ(this->count_cells(tiles), 4845);
    EXPECT_EQ(tiles.size(), 48);
}

TYPED_TEST(TileMappingTest, asymmetricImpl1)
{
    using TileSet_t     = TestFixture::TileSet_t;
    using TileSetInfo_t = TestFixture::info_t;
    auto const& [box, shape, tile_size]
        = TileSetInfo_t{{{0, 0, 0}, {18, 16, 14}}, {3, 3, 3}, {5, 5, 5}};

    auto const tiles = Tiler<TileSet_t, 1>{{box, shape, tile_size}}.f();

    EXPECT_TRUE(this->all_in(tiles, box));
    EXPECT_FALSE(any_overlaps_in(tiles, [](auto const& tile) { return tile; }));
    EXPECT_EQ(this->count_cells(tiles), product(box.shape()));
    EXPECT_EQ(this->count_cells(tiles), 4845);
    EXPECT_EQ(tiles.size(), product(shape));
    EXPECT_EQ(tiles.size(), 27);
}

TYPED_TEST(TileMappingTest, asymmetricImpl1_2)
{
    using TileSet_t     = TestFixture::TileSet_t;
    using TileSetInfo_t = TestFixture::info_t;

    // not used in impl 1 mapper
    std::array<std::uint32_t, TileSet_t::dimension> const tiles_per_dim{0, 0, 0};

    for (std::uint8_t i = 10; i < 11; ++i)
        for (std::uint8_t j = 10; j < 11; ++j)
            for (std::uint8_t k = 10; k < 11; ++k)
            {
                auto const& [box, shape, tile_size]
                    = TileSetInfo_t{{{0, 0, 0}, {i, j, k}}, tiles_per_dim, {4, 4, 4}};
                auto const tiles = Tiler<TileSet_t, 1>{{box, shape, tile_size}}.f();

                EXPECT_TRUE(this->all_in(tiles, box));
                EXPECT_FALSE(any_overlaps_in(tiles, [](auto const& tile) { return tile; }));
                EXPECT_EQ(this->count_cells(tiles), product(box.shape()));
                this->print(tiles);
                break;
            }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
