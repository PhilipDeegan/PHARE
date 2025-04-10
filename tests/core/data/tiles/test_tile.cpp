
#include "test_tile.hpp"

using DimTiles = testing::Types<TileSet<TileMock<1>>, TileSet<TileMock<2>>, TileSet<TileMock<3>>>;

TYPED_TEST_SUITE(TileTestBoxShapeNotMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTestBoxShapeMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTest, DimTiles);


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, expectedNbrOfTilesPerDimToCoverTheBox)
{
    auto const& shape = this->tileSet.shape();
    for (auto i = 0u; i < this->dimension; ++i)
        EXPECT_EQ(shape[i], 14);
}

TYPED_TEST(TileTestBoxShapeMultipleTileSize, cluserSetSizeIsCorrect)
{
    EXPECT_EQ(this->tileSet.size(), std::pow(12, this->dimension));
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, sumOfTileSurfacesEqualsTileSetBoxSurface)
{
    auto surface = std::accumulate(std::begin(this->tileSet), std::end(this->tileSet), 0u,
                                   [&](auto acc, auto const& tile) { return acc + tile.size(); });
    EXPECT_EQ(surface, this->box.size());
}



TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, tileHasNoOverlapWithOthers)
{
    auto constexpr dim = TypeParam::dimension;
    for (auto const& tile : this->tileSet)
    {
        for (auto const& other : this->tileSet)
        {
            if (&tile != &other)
            {
                auto overlap = tile * other;
                EXPECT_FALSE(overlap.has_value());
            }
        }
    }
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, retrieveTilesFromBoxOverlap)
{
    auto constexpr dim = TypeParam::dimension;
    Box<int, dim> selection_box{ConstArray<int, dim>(11), ConstArray<int, dim>(34)};

    auto expected_nbr = std::pow(7, this->dimension);
    auto overlapeds   = this->tileSet.overlaped_with(selection_box);
    EXPECT_EQ(overlapeds.size(), expected_nbr);

    auto completes   = 0.;
    auto incompletes = 0.;
    for (auto const& overlaped : overlapeds)
    {
        auto const& [is_complete, tile] = overlaped;
        if (is_complete)
            ++completes;
        else
            ++incompletes;
    }
    EXPECT_EQ(completes, std::pow(5, dim));
    EXPECT_EQ(incompletes, std::pow(7, dim) - std::pow(5, dim));
}


TYPED_TEST(TileTest, cannotCreateTileWithTileSizeBiggerThanBox)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(5)};
    auto const tile_size = PHARE::core::ConstArray<std::size_t, dim>(7); // larger than box shape
    // EXPECT_THROW(std::make_unique<TypeParam>(box, tile_size), std::runtime_error);
}


TYPED_TEST(TileTest, linksAreValid)
{
    if constexpr (TypeParam::dimension == 3)
    {
        TypeParam::build_links(this->tileSet);
        std::unordered_set<TileMock<3>*> links;
        std::function<void(TileMock<3>*)> traverse = [&](TileMock<3>* tile) {
            if (!tile || links.count(tile))
                return;
            links.emplace(tile);
            for (auto const link : tile->links())
                traverse(link);
        };
        traverse(&this->tileSet[0]);
        EXPECT_EQ(std::pow(12, 3), links.size());
    }
}

TYPED_TEST(TileTest, canTraverseLinks)
{
    if constexpr (TypeParam::dimension == 3)
    {
        TypeParam::build_links(this->tileSet);
        traverse_ghost_boundary_tiles(this->tileSet, [](auto& tile) {
            // do something
        });
    }
}


TYPED_TEST(TileTest, canTraverseNeighbours)
{
    if constexpr (TypeParam::dimension == 3)
    {
        TypeParam::build_links(this->tileSet);
        traverse_tile_neighbours(this->tileSet, *this->tileSet[0].link(6), [](auto& tile) {
            // do something
        });
    }
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, canRetrieveTileFromCell)
{
    auto constexpr dim = TypeParam::dimension;
    auto tile          = [&]() {
        if constexpr (dim == 1)
            return this->tileSet.at(13);
        else if constexpr (dim == 2)
            return this->tileSet.at(13, 13);
        else if constexpr (dim == 3)
            return this->tileSet.at(13, 13, 13);
    }();
    auto const expected_box = Box<int, dim>{ConstArray<int, dim>(12), ConstArray<int, dim>(15)};
    EXPECT_TRUE(*tile == expected_box);
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, useTileSetView)
{
    auto view  = this->tileSet.make_view();
    auto shape = view.shape();
    for (auto const& tile : view)
    {
        EXPECT_LE(tile.size(), std::pow(4, this->dimension));
    }
}


template<std::size_t dim>
struct PatchDataMock
{
    using BoxND     = Box<int, dim>;
    using TileSetND = TileSet<BoxND>;
    Box<int, dim> box;
    TileSetND tileSet;

    PatchDataMock(BoxND box_,
                  std::array<std::size_t, dim> const& tile_size = ConstArray<std::size_t, dim>(4))
        : box{box_}
        , tileSet{box, tile_size}
    {
    }
};

template<std::size_t dim>
struct RankDataMock
{
    RankDataMock()
    {
        if constexpr (dim == 1)
        {
            patchdatas.emplace_back(Box<int, dim>{{0}, {24}});
            patchdatas.emplace_back(Box<int, dim>{{33}, {65}});
            patchdatas.emplace_back(Box<int, dim>{{102}, {134}});
        }
        else if constexpr (dim == 2)
        {
            patchdatas.emplace_back(Box<int, dim>{{0, 10}, {24, 58}});
            patchdatas.emplace_back(Box<int, dim>{{33, 57}, {65, 83}});
            patchdatas.emplace_back(Box<int, dim>{{102, 99}, {134, 128}});
        }
        else if constexpr (dim == 3)
        {
            patchdatas.emplace_back(Box<int, dim>{{0, 10, 20}, {24, 58, 78}});
            patchdatas.emplace_back(Box<int, dim>{{33, 57, 77}, {65, 83, 98}});
            patchdatas.emplace_back(Box<int, dim>{{102, 99, 120}, {134, 128, 138}});
        }
    }
    std::vector<PatchDataMock<dim>> patchdatas;
};


TEST(TileSetViewSpan, fromManyPatches)
{
    RankDataMock<2> rankData;
    std::vector<TileSetView<Box<int, 2>>> views;

    // now making views
    for (auto /*const? make_view fails if so...*/& patch : rankData.patchdatas)
    {
        views.push_back(patch.tileSet.make_view());
    }
}


TEST(rankData, canCreateRankData)
{
    RankDataMock<1> rankData;
    RankDataMock<2> rankData2;
    RankDataMock<3> rankData3;
}




TYPED_TEST(TileTestBoxShapeMultipleTileSize, InnerTileHaveCorrectNbrOfNeighbors)
{
    // build the tileSet for a given box
    auto constexpr tile_size = 4u;
    auto constexpr dim       = TypeParam::dimension;
    using BoxND              = Box<int, dim>;

    auto const lower = ConstArray<int, dim>(0);
    auto const upper = ConstArray<int, dim>(51);
    BoxND box{lower, upper};
    TileSet<BoxND> tileSet{box, ConstArray<std::size_t, dim>(tile_size)};

    auto const inner_tiles = this->tileSet.inner_tiles();

    auto constexpr expected_neighbor_nbr = [&]() {
        if constexpr (dim == 1)
            return 2;
        else if constexpr (dim == 2)
            return 8;
        else if constexpr (dim == 3)
            return 26;
    }();

    for (auto const& tile : inner_tiles)
    {
        // looping over all cells of the tileSet box +1 ghost cell
        // we should see 9 different pointers, that is (2,8,26) neighbors + the tile itself
        // we discard the tile itself
        std::unordered_set<BoxND*> neighbors;
        BoxND ghost_box{*tile};
        ghost_box.grow(1);
        for (auto const& cell : ghost_box)
        {
            auto tile_ptr = [&]() {
                if constexpr (dim == 1)
                    return this->tileSet.at(cell[0]);
                else if constexpr (dim == 2)
                    return this->tileSet.at(cell[0], cell[1]);
                else
                    return this->tileSet.at(cell[0], cell[1], cell[2]);
            }();
            if (!isIn(cell, static_cast<BoxND>(*tile)))
                neighbors.insert(tile_ptr);
        }
        EXPECT_EQ(neighbors.size(), expected_neighbor_nbr);
    }
}


TYPED_TEST(TileTestBoxShapeMultipleTileSize, gettingBorderTiles)
{
    auto border_tiles = this->tileSet.border_tiles();
    auto inner_tiles  = this->tileSet.inner_tiles();
    auto total_tiles  = this->tileSet.size();
    EXPECT_EQ(border_tiles.size(), total_tiles - inner_tiles.size());
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
