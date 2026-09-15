
#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_converter.hpp"
#include "core/data/particles/particle_array_comparator.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"


namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{14});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{1000});


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
};


template<typename Param>
struct ParticleArrayConstructionTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static sim_opts
        = SimOpts{.dimension = dim, .layout_mode = layout_mode, .alloc_mode = alloc_mode};

    using GridLayout_t = TestGridLayout<typename PHARE_Types<sim_opts>::Hybrid::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    GridLayout_t layout{cells};

    ParticleArray_t setup_particles() const // test movable
    {
        auto ps = make_particles<ParticleArray_t>(layout);
        add_particles(ps, layout.AMRBox(), ppc);
        delta_disperse(ps);
        vary_velocity(ps, -6, 6);
        return ps;
    }
};



// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

    TestParam<1, LayoutMode::AoS, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>

>;
// clang-format on



TYPED_TEST_SUITE(ParticleArrayConstructionTest, Permutations_t, );


// 2^dim "corner" neighbours only (every coordinate moves by +/-1), ordered so that
// treating -1 as bit 0 and +1 as bit 1 counts up from (-1,-1,-1) to (1,1,1).
template<std::size_t dim>
auto corner_offsets()
{
    std::array<std::array<int, dim>, (std::size_t{1} << dim)> offsets{};
    for (std::size_t idx = 0; idx < offsets.size(); ++idx)
        for (std::size_t d = 0; d < dim; ++d)
            offsets[idx][d] = ((idx >> (dim - 1 - d)) & 1) ? 1 : -1;
    return offsets;
}

template<typename ICell>
auto add_icell(ICell const& a, ICell const& b)
{
    ICell out;
    for (std::size_t d = 0; d < a.size(); ++d)
        out[d] = a[d] + b[d];
    return out;
}


// Function templates can't be partially specialised (fixing layout_mode while leaving
// Particles deduced), so the per-layout dispatch lives on a class template instead,
// which can be. One specialisation per layout_mode, each still generic over Particles.
template<auto layout_mode>
struct MoveParticles;

template<>
struct MoveParticles<LayoutMode::AoS>
{
    // no registration needed, flat array has no cell/tile structure to maintain
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        std::size_t counter = 0;
        for (auto& p : particles)
        {
            p.iCell() = add_icell(p.iCell(), offsets[counter % offsets.size()]);
            ++counter;
        }
    }
};

template<>
struct MoveParticles<LayoutMode::AoSMapped>
{
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        std::size_t counter = 0;
        for (std::size_t idx = 0; idx < particles.size(); ++idx)
        {
            auto const newcell
                = add_icell(particles[idx].iCell(), offsets[counter % offsets.size()]);
            ++counter;
            particles.change_icell(newcell, idx);
        }
    }
};

template<>
struct MoveParticles<LayoutMode::AoSPC>
{
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        for (auto const& bix : particles.local_box())
        {
            auto& cell_particles = particles(bix);
            auto const n         = cell_particles.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto& p            = cell_particles[i];
                auto const oldcell = p.iCell();
                auto const newcell = add_icell(oldcell, offsets[counter % offsets.size()]);
                ++counter;
                p.iCell() = newcell;
                ParticleTracker<dim> const pt{oldcell};
                particles.template move_check<ParticleType::Domain>(pt, i, p);
            }
        }
        particles.template on_moved<ParticleType::Domain>();
    }
};


template<>
struct MoveParticles<LayoutMode::AoSPCTS>
{
    // registration goes through the span (move_check lives on PCTileSetSpan); the view's
    // gap arrays alias the vector's, so on_moved() on the vector sees every registration
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        auto view           = particles.view();
        for (auto& tile : view())
        {
            // tile_cell names the tile PHYSICALLY holding these particles (this same
            // loop iteration), not view.local_tile_cell(oldcell)'s clamp-based
            // "owner" — see MoveLevelGhostParticles<AoSPCTS> and
            // multi_boris.hpp::per_tile for why these must match
            auto const tile_cell = view.local_cell(tile.lower);
            auto& tile_particles = tile();
            for (auto const& bix : tile_particles.local_box())
            {
                auto& cell_particles = tile_particles(bix);
                auto const n         = cell_particles.size();
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto& p            = cell_particles[i];
                    auto const oldcell = p.iCell();
                    auto const newcell = add_icell(oldcell, offsets[counter % offsets.size()]);
                    ++counter;
                    p.iCell() = newcell;
                    TiledParticleTracker<dim> const pt{oldcell, tile_cell};
                    view.template move_check<ParticleType::Domain>(pt, i, p);
                }
            }
        }
        particles.template on_moved<ParticleType::Domain>();
    }
};


template<typename Particles>
void move_particles(Particles& particles)
{
    MoveParticles<Particles::layout_mode>::apply(particles, corner_offsets<Particles::dimension>());
}


// level ghost arrays: particles live in the ghost layer; movers re-bucket anywhere
// inside the ghost box (including into the domain), ghost-box leavers are deleted at
// sync. only per-cell layouts register level ghost moves — other layouts skip the test
template<auto layout_mode>
struct MoveLevelGhostParticles;

template<>
struct MoveLevelGhostParticles<LayoutMode::AoS>
{
    // no registration: move directly, then apply the deletion contract by hand
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const& ghost_box)
    {
        std::size_t counter = 0;
        for (auto& p : particles)
        {
            p.iCell() = add_icell(p.iCell(), offsets[counter % offsets.size()]);
            ++counter;
        }
        std::erase_if(particles.vector(),
                      [&](auto const& p) { return not isIn(p.iCell(), ghost_box); });
    }
};

template<>
struct MoveLevelGhostParticles<LayoutMode::AoSPC>
{
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const&)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        for (auto const& bix : particles.local_box())
        {
            auto& cell_particles = particles(bix);
            auto const n         = cell_particles.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto& p            = cell_particles[i];
                auto const oldcell = p.iCell();
                p.iCell()          = add_icell(oldcell, offsets[counter % offsets.size()]);
                ++counter;
                ParticleTracker<dim> const pt{oldcell};
                particles.template move_check<ParticleType::LevelGhost>(pt, i, p);
            }
        }
        particles.template on_moved<ParticleType::LevelGhost>();
    }
};

template<>
struct MoveLevelGhostParticles<LayoutMode::AoSPCTS>
{
    // registration goes through the span, as for the domain mover above
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const&)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        auto view           = particles.view();
        for (auto& tile : view())
        {
            // tile_cell names the tile PHYSICALLY holding these particles (this same
            // loop iteration), not view.local_tile_cell(oldcell)'s clamp-based
            // "owner" tile — adjacent tiles' grown ghost boxes legitimately overlap at
            // shared tile-to-tile ghost cells (particle_array_refiner.hpp writes a
            // duplicate physical copy into every overlapping tile there), so the
            // clamp owner of a given cell can differ from the tile that physically
            // holds a given particle. move_check's LevelGhost branch relies on
            // pt.tile_cell reflecting physical residency (see its comments).
            auto const tile_cell = view.local_cell(tile.lower);
            auto& tile_particles = tile();
            for (auto const& bix : tile_particles.local_box())
            {
                auto& cell_particles = tile_particles(bix);
                auto const n         = cell_particles.size();
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto& p            = cell_particles[i];
                    auto const oldcell = p.iCell();
                    p.iCell()          = add_icell(oldcell, offsets[counter % offsets.size()]);
                    ++counter;
                    auto const pt
                        = make_particle_tracker<LayoutMode::AoSPCTS, ParticleType::LevelGhost, dim>(
                            oldcell, tile_cell);
                    view.template move_check<ParticleType::LevelGhost>(pt, i, p);
                }
            }
        }
        particles.template on_moved<ParticleType::LevelGhost>();
    }
};


// every per-cell bucket must only hold particles whose iCell maps to that bucket
template<typename Particles>
void check_cell_buckets(Particles const& particles)
{
    using enum LayoutMode;
    auto constexpr static dim = Particles::dimension;

    if constexpr (any_in(Particles::layout_mode, AoSPC, AoSPCTS))
    {
        auto const check = [](auto const& cps, std::array<std::uint32_t, dim> const& cell) {
            for (auto const& p : cps(cell))
                EXPECT_TRUE(array_equals(cps.local_cell(p.iCell()), cell))
                    << "particle in wrong cell bucket: " << Point{p.iCell()};
        };

        if constexpr (Particles::layout_mode == AoSPC)
        {
            for (auto const& bix : particles.local_box())
                check(particles, bix);
        }
        else
        {
            for (auto const& tile : particles())
                for (auto const& bix : tile().local_box())
                    check(tile(), bix);
        }
    }
}


// only particles outside the patch domain box may sit outside their tile's box —
// a domain particle must never be stored in a tile's ghost cells. LevelGhost is
// exempt: a level ghost particle may legitimately stay self-contained in the tile
// whose own ghost box it drifted within, even when its iCell() lands in a domain cell
// that "belongs" to a neighbouring tile — move_in_domain's later scan is tile-
// ownership-agnostic, it only needs the particle discoverable in some tile's own
// bucket at the right AMR cell (see particle_array_pc_ts.hpp move_check comments)
template<auto particle_type = ParticleType::Domain, typename Particles>
void check_tile_ownership(Particles const& particles)
{
    using enum LayoutMode;

    if constexpr (particle_type == ParticleType::LevelGhost)
        return;
    else if constexpr (any_in(Particles::layout_mode, AoSPCTS))
    {
        auto const check = [&](auto const& tile, auto const& p) {
            if (isIn(p.iCell(), particles.box()))
            {
                EXPECT_TRUE(isIn(p.iCell(), tile))
                    << "domain particle in tile ghost cells: " << Point{p.iCell()}
                    << " tile: " << tile;
            }
        };

        for (auto const& tile : particles())
        {
            auto const& cps = tile();
            for (auto const& bix : cps.local_box())
                for (auto const& p : cps(bix))
                    check(tile, p);
        }
    }
}



TYPED_TEST(ParticleArrayConstructionTest, test_move_sync_works)
{
    using ParticleArray_t    = TestFixture::ParticleArray_t;
    using AoSParticleArray_t = AoSParticleArray<TestFixture::dim>;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles(particles, this->layout.AMRBox(), ppc);

    // independent AoS ground truth, same initial particles, moved directly (no
    // registration needed for AoS) so it never depends on the layout under test
    auto reference = convert_particles<AoSParticleArray_t>(particles, this->layout);

    move_particles(particles);
    move_particles(reference);

    check_tile_ownership(particles);
    check_cell_buckets(particles);

    // particles have moved out of AMRBox into the first ghost layer; sorting must use the
    // grown box or distinct out-of-box cells alias to the same flat index and, with all
    // deltas identical, the two arrays end up in different orders
    auto const sort_box = grow(this->layout.AMRBox(), 1);
    auto converted      = convert_particles<AoSParticleArray_t>(particles, this->layout);
    sort_particles(converted, sort_box);
    sort_particles(reference, sort_box);

    auto const report = compare_particles(reference, converted);
    EXPECT_TRUE(report) << report.why();
}


// level ghost arrays fill the ghost layer; after a move they may re-bucket anywhere in
// the ghost box (including into the domain) and ghost-box leavers must be deleted
TYPED_TEST(ParticleArrayConstructionTest, test_level_ghost_move_sync_works)
{
    using ParticleArray_t    = TestFixture::ParticleArray_t;
    using AoSParticleArray_t = AoSParticleArray<TestFixture::dim>;
    using enum LayoutMode;

    auto constexpr static dim = TestFixture::dim;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    if constexpr (ParticleArray_t::alloc_mode != AllocatorMode::CPU
                  or not any_in(ParticleArray_t::layout_mode, AoS, AoSPC, AoSPCTS))
        GTEST_SKIP() << "level ghost move_check unsupported for this layout";
    else
    {
        auto constexpr static ghosts = TestFixture::GridLayout_t::options.particle_ghost_width;
        auto const ghost_box         = grow(this->layout.AMRBox(), ghosts);

        auto particles = make_particles<ParticleArray_t>(this->layout);
        add_ghost_particles(particles, this->layout.AMRBox(), ppc, ghosts);

        auto reference = convert_particles<AoSParticleArray_t>(particles, this->layout);

        auto const offsets = corner_offsets<dim>();
        MoveLevelGhostParticles<ParticleArray_t::layout_mode>::apply(particles, offsets, ghost_box);
        MoveLevelGhostParticles<LayoutMode::AoS>::apply(reference, offsets, ghost_box);

        check_tile_ownership<ParticleType::LevelGhost>(particles);
        check_cell_buckets(particles);

        auto converted = convert_particles<AoSParticleArray_t>(particles, this->layout);
        sort_particles(converted, ghost_box);
        sort_particles(reference, ghost_box);

        EXPECT_EQ(reference.size(), converted.size());
        auto const report = compare_particles(reference, converted);
        EXPECT_TRUE(report) << report.why();
    }
}


// Reproduces the AoSPCTS cross-tile LevelGhost bug (see mkn/now notes): the L0->L1
// coarse-to-fine refiner (particle_array_refiner.hpp:328-340) deliberately writes a
// duplicate physical copy of a level ghost particle into EVERY tile whose own grown
// ghost box reaches a shared tile-to-tile ghost cell, not just the single tile
// TileSet::at() clamp-designates as "owner" of that cell (used for field-side
// reduce_into_, which already expects and sums duplicate contributions). But
// move_check/sync_moved/copy_in (particle_array_pc_ts.hpp) all key their cross-tile
// gap/capacity bookkeeping by a single shared per-AMR-cell array (add_into_/gap_idx_/
// gaps_/cell_size_), implicitly assuming exactly one physical bucket per cell — the
// clamp-owner's. A second physical bucket in a non-owning tile silently corrupts that
// shared bookkeeping once both copies move, since move_check has no way to tell which
// physical bucket a registered gap index belongs to (see move_check's fallthrough:
// `(*particles_.at(pt.tile_cell))().move_check(pt, idx, particle)` always resolves to
// the clamp owner, never the tile that physically holds the data).
//
// This test manufactures that duplication directly (without needing the full
// refiner/hierarchy machinery) by writing a second copy of an existing level ghost
// particle straight into a non-owning neighbour tile's own per-cell storage, then runs
// the normal move + on_moved<LevelGhost>() cycle used elsewhere in this file.
TYPED_TEST(ParticleArrayConstructionTest, test_level_ghost_cross_tile_duplicate)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;
    using enum LayoutMode;

    auto constexpr static dim = TestFixture::dim;

    // cross-tile ghost-cell overlap needs a wall (one dim) crossed with the level
    // ghost layer (another dim) - can't happen in 1D, and only AoSPCTS duplicates
    if constexpr (ParticleArray_t::layout_mode != AoSPCTS or dim < 2)
        GTEST_SKIP() << "cross-tile level ghost duplication only applies to multi-dim AoSPCTS";
    else
    {
        auto constexpr static ghosts = TestFixture::GridLayout_t::options.particle_ghost_width;
        auto const& amr_box          = this->layout.AMRBox();
        auto const ghost_box         = grow(amr_box, ghosts);

        auto particles = make_particles<ParticleArray_t>(this->layout);
        add_ghost_particles(particles, amr_box, ppc, ghosts);

        // find a level ghost cell (outside the whole patch) reachable by more than
        // one tile's own grown ghost box - a shared tile-to-tile ghost cell the
        // refiner would legitimately write more than one physical copy into
        using Tile_t       = typename ParticleArray_t::Tile_t;
        Tile_t* owner_tile = nullptr;
        Tile_t* dup_tile   = nullptr;
        std::array<int, dim> dup_cell{};

        for (auto& tile : particles())
        {
            for (auto const& cell : grow(tile, ghosts))
            {
                if (isIn(cell, amr_box))
                    continue; // only true level ghost cells (outside the whole patch)

                auto* owner = particles().at(cell);
                if (owner != &tile)
                {
                    owner_tile = owner;
                    dup_tile   = &tile;
                    dup_cell   = cell.toArray();
                    break;
                }
            }
            if (dup_tile)
                break;
        }

        ASSERT_NE(dup_tile, nullptr)
            << "expected at least one tile-to-tile shared level ghost cell for this tiling";

        // duplicate one already-existing particle at that cell directly into the
        // non-owning tile's own per-cell storage, exactly as the refiner does
        auto const owner_lcell = (*owner_tile)().local_cell(dup_cell);
        ASSERT_GT((*owner_tile)()(owner_lcell).size(), std::size_t{0});
        auto const p = (*owner_tile)()(owner_lcell)[0];

        auto& dup_pc         = (*dup_tile)();
        auto const dup_lcell = dup_pc.local_cell(dup_cell);
        dup_pc(dup_lcell).reserve(dup_pc(dup_lcell).size() + 1);
        dup_pc(dup_lcell).emplace_back(p);

        particles.template on_appended<ParticleType::LevelGhost>();

        auto const offsets = corner_offsets<dim>();
        MoveLevelGhostParticles<AoSPCTS>::apply(particles, offsets, ghost_box);

        check_tile_ownership<ParticleType::LevelGhost>(particles);
        check_cell_buckets(particles);
    }
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
