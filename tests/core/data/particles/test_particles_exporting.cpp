
// exporters (delete_particles_(not_)in, move_in_domain, move_in_ghost_layer) per layout and
//  allocator, against a reference computed from each operation's definition on a plain vector

#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_exporter.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_selector.hpp"

#include "tests/core/data/particles/test_particles.hpp"

#include "gtest/gtest.h"

#include <tuple>
#include <vector>
#include <algorithm>


namespace PHARE::core
{

template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;

    using ParticleArray_t = ParticleArray<ParticleArrayOptions{dim, lm, StorageMode::VECTOR, am}>;
};


template<typename Param>
struct ParticlesExportingTest : public ::testing::Test
{
    auto constexpr static dim = Param::dim;
    using ParticleArray_t     = Param::ParticleArray_t;
    using Particle_t          = Particle<dim>;
    using Box_t               = Box<int, dim>;
    using Ref_t               = std::vector<Particle_t>;

    std::size_t constexpr static ghost_cells = 2;
    std::size_t constexpr static ppc         = 3;
    int constexpr static cells               = 12; // more than one tile per dimension

    Box_t const domain{ConstArray<int, dim>(0), ConstArray<int, dim>(cells - 1)};
    Box_t const ghost_box = grow(domain, ghost_cells);

    // layouts whose particles can have their cell changed in place, without a remap: as after a
    //  push, before the exporters run. Tiled ones then hold particles a cell past their ghost box
    bool constexpr static can_push
        = any_in(Param::layout_mode, LayoutMode::AoS, LayoutMode::AoSTS);

    // AoSMapped declares the move exporters but has no definitions
    bool constexpr static can_move = Param::layout_mode != LayoutMode::AoSMapped;

    // unique per id, so moved / deleted particles are told apart
    static Particle_t make(std::array<int, dim> const& icell, std::size_t const id)
    {
        auto const frac = [](double const x) { return x - static_cast<int>(x); };
        std::array<double, dim> delta;
        for (std::size_t d = 0; d < dim; ++d)
            delta[d] = frac(id * .6180339887 + d * .3247179572);
        return Particle_t{.001 + id * 1e-9, icell, delta, {id * 1e-3, id * 2e-3, id * 3e-3}};
    }

    // the in place "push": by delta, so the same particles move in every array
    static void push(auto& p)
    {
        if (p.delta()[0] < .2)
            p.iCell()[0] -= 1;
        else if (p.delta()[0] > .8)
            p.iCell()[0] += 1;
    }

    auto make_array() const { return make_particles<ParticleArray_t>(domain, ghost_cells); }

    // the array under test and its reference, ppc particles per ghost box cell
    auto fill(bool const pushed) const
    {
        auto ps = make_array();
        Ref_t ref;
        core::add_particles<ParticleType::Domain>(ps, ghost_box, ppc, make);
        std::size_t id = 0;
        for_each_point(ghost_box, [&](auto const& bix) {
            for (std::size_t i = 0; i < ppc; ++i)
                ref.emplace_back(make(*bix, id++));
        });
        if (pushed)
        {
            per_particle(ps, [](auto& p) { push(p); });
            for (auto& p : ref)
                push(p);
        }
        EXPECT_EQ(ps.size(), ref.size());
        return std::make_tuple(std::move(ps), std::move(ref));
    }

    static auto key(auto const& p)
    {
        return std::make_tuple(p.iCell(), p.delta(), p.v(), p.weight());
    }

    static auto sorted(Ref_t ref)
    {
        std::sort(ref.begin(), ref.end(),
                  [](auto const& a, auto const& b) { return key(a) < key(b); });
        return ref;
    }

    template<typename Array_t>
    static auto gather(Array_t& ps)
    {
        Ref_t out;
        per_particle(ps, [&](auto const& p) {
            out.emplace_back(p.weight(), p.iCell(), p.delta(), p.v());
        });
        return sorted(std::move(out));
    }

    template<typename Array_t>
    static void expect_same(Array_t& ps, Ref_t const& ref, std::string const& what)
    {
        EXPECT_EQ(ps.size(), ref.size()) << what;
        auto const got = gather(ps);
        auto const exp = sorted(ref);
        ASSERT_EQ(got.size(), exp.size()) << what;
        for (std::size_t i = 0; i < got.size(); ++i)
        {
            ASSERT_EQ(key(got[i]), key(exp[i])) << what << " at " << i;
        }
        if constexpr (requires { ps.is_consistent(); })
        {
            EXPECT_TRUE(ps.is_consistent()) << what;
        }
    }

    // runs check for the settled fill, and for the pushed one where the layout allows it
    void for_fills(auto&& check) const
    {
        check(false);
        if constexpr (can_push)
            check(true);
    }

    static auto keep_if(Ref_t const& ref, auto&& pred)
    {
        Ref_t out;
        std::copy_if(ref.begin(), ref.end(), std::back_inserter(out), pred);
        return out;
    }

    static bool in_any(auto const& p, auto const& boxes)
    {
        return std::any_of(boxes.begin(), boxes.end(),
                           [&](auto const& box) { return isIn(p, box); });
    }

    // splits box in two along the first dimension
    static auto halves(Box_t const& box)
    {
        auto lo = box, hi = box;
        auto const mid = (box.lower[0] + box.upper[0]) / 2;
        lo.upper[0]    = mid;
        hi.lower[0]    = mid + 1;
        return std::vector<Box_t>{lo, hi};
    }

    static std::string label(bool const pushed) { return pushed ? "pushed" : "settled"; }
};


// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
    TestParam<1, LayoutMode::AoS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSCMTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSCMTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSCMTS, AllocatorMode::CPU>

PHARE_WITH_GPU(
   ,TestParam<1, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED>
)

>;
// clang-format on

TYPED_TEST_SUITE(ParticlesExportingTest, Permutations_t, );


TYPED_TEST(ParticlesExportingTest, delete_particles_not_in_box)
{
    this->for_fills([&](bool const pushed) {
        auto [ps, ref] = this->fill(pushed);
        delete_particles_not_in(ps, this->domain);
        this->expect_same(ps, this->keep_if(ref, [&](auto const& p) {
            return isIn(p, this->domain);
        }), this->label(pushed));
    });
}

TYPED_TEST(ParticlesExportingTest, delete_particles_not_in_boxes)
{
    this->for_fills([&](bool const pushed) {
        auto [ps, ref]    = this->fill(pushed);
        auto const boxes  = this->halves(this->domain);
        delete_particles_not_in(ps, boxes);
        this->expect_same(ps, this->keep_if(ref, [&](auto const& p) {
            return this->in_any(p, boxes);
        }), this->label(pushed));
    });
}

// the ghost box as the box: border tiles grown by the ghost width fit inside it, but pushed
//  particles sit a cell past it
TYPED_TEST(ParticlesExportingTest, delete_particles_not_in_ghost_box)
{
    this->for_fills([&](bool const pushed) {
        auto [ps, ref] = this->fill(pushed);
        delete_particles_not_in(ps, this->ghost_box);
        this->expect_same(ps, this->keep_if(ref, [&](auto const& p) {
            return isIn(p, this->ghost_box);
        }), this->label(pushed));
    });
}

TYPED_TEST(ParticlesExportingTest, delete_particles_not_in_ghost_boxes)
{
    this->for_fills([&](bool const pushed) {
        auto [ps, ref]   = this->fill(pushed);
        auto const boxes = this->halves(this->ghost_box);
        delete_particles_not_in(ps, boxes);
        this->expect_same(ps, this->keep_if(ref, [&](auto const& p) {
            return this->in_any(p, boxes);
        }), this->label(pushed));
    });
}

TYPED_TEST(ParticlesExportingTest, delete_particles_in_boxes)
{
    using Box_t = TestFixture::Box_t;
    auto const& domain = this->domain;
    // one box inside the domain (as a non owned box), one across the upper ghost layer
    std::vector<Box_t> const boxes{
        Box_t{domain.lower + 2, domain.lower + 4},
        Box_t{domain.upper - 1, domain.upper + 1},
    };

    this->for_fills([&](bool const pushed) {
        auto [ps, ref] = this->fill(pushed);
        delete_particles_in(ps, boxes);
        this->expect_same(ps, this->keep_if(ref, [&](auto const& p) {
            return !this->in_any(p, boxes);
        }), this->label(pushed));
    });
}

TYPED_TEST(ParticlesExportingTest, move_in_domain)
{
    if constexpr (!TestFixture::can_move)
        GTEST_SKIP() << "no move exporters for this layout";
    else
    {
    this->for_fills([&](bool const pushed) {
        auto [src, ref] = this->fill(pushed);
        auto dst        = this->make_array();
        move_in_domain(dst, src, this->domain);

        auto const moves = [&](auto const& p) { return isIn(p, this->domain); };
        this->expect_same(src, this->keep_if(ref, [&](auto const& p) { return !moves(p); }),
                          "src " + this->label(pushed));
        this->expect_same(dst, this->keep_if(ref, moves), "dst " + this->label(pushed));
    });
    }
}

TYPED_TEST(ParticlesExportingTest, move_in_ghost_layer_box)
{
    if constexpr (!TestFixture::can_move)
        GTEST_SKIP() << "no move exporters for this layout";
    else
    {
    this->for_fills([&](bool const pushed) {
        auto [src, ref] = this->fill(pushed);
        auto dst        = this->make_array();
        move_in_ghost_layer(dst, src, this->domain, this->ghost_box);

        auto const moves = [&](auto const& p) {
            return isIn(p, this->ghost_box) and !isIn(p, this->domain);
        };
        this->expect_same(src, this->keep_if(ref, [&](auto const& p) { return !moves(p); }),
                          "src " + this->label(pushed));
        this->expect_same(dst, this->keep_if(ref, moves), "dst " + this->label(pushed));
    });
    }
}

TYPED_TEST(ParticlesExportingTest, move_in_ghost_layer_boxes)
{
    if constexpr (!TestFixture::can_move)
        GTEST_SKIP() << "no move exporters for this layout";
    else
    {
    auto const ghost_boxes = this->halves(this->ghost_box);

    this->for_fills([&](bool const pushed) {
        auto [src, ref] = this->fill(pushed);
        auto dst        = this->make_array();
        move_in_ghost_layer(dst, src, this->domain, ghost_boxes);

        auto const moves = [&](auto const& p) {
            return this->in_any(p, ghost_boxes) and !isIn(p, this->domain);
        };
        this->expect_same(src, this->keep_if(ref, [&](auto const& p) { return !moves(p); }),
                          "src " + this->label(pushed));
        this->expect_same(dst, this->keep_if(ref, moves), "dst " + this->label(pushed));
    });
    }
}


// appending from a host AoSMapped array (as the refiner does), into the tile owning each cell
TYPED_TEST(ParticlesExportingTest, append_particles_from_aos_mapped)
{
    using Src_t = ParticleArray<
        ParticleArrayOptions{TestFixture::dim, LayoutMode::AoSMapped, StorageMode::VECTOR}>;

    if constexpr (!TestFixture::can_move) // no AoSMapped to AoSMapped appender either
        GTEST_SKIP() << "no appender for this layout";
    else
    {
        auto src = make_particles<Src_t>(this->domain, this->ghost_cells);
        core::add_particles<ParticleType::Domain>(src, this->ghost_box, this->ppc,
                                                  TestFixture::make);
        auto [dst, ref] = this->fill(false);
        append_particles<ParticleType::Domain>(src, dst);

        auto expected = ref;
        expected.insert(expected.end(), ref.begin(), ref.end()); // the same particles twice
        this->expect_same(dst, expected, "append");
    }
}


// selections as ParticlesData does them: into a flat array (the pack buffer), shifted (periodic),
//  into the same layout (patch copy), and counts. Domain selections only see a tile's own box,
//  so their boxes stay inside the domain

TYPED_TEST(ParticlesExportingTest, select_particles_into_flat)
{
    using Flat_t = ParticleArray<ParticleArrayOptions{
        TestFixture::dim, LayoutMode::AoS, StorageMode::VECTOR, TypeParam::alloc_mode}>;
    using Box_t = TestFixture::Box_t;

    auto [src, ref] = this->fill(false);
    Box_t const box{this->domain.lower + 1, this->domain.upper - 2};
    Flat_t dst;
    select_particles<ParticleType::Domain>(src, dst, box);
    this->expect_same(dst, this->keep_if(ref, [&](auto const& p) { return isIn(p, box); }),
                      "select");
    EXPECT_EQ(count_particles<ParticleType::Domain>(src, box), dst.size());
}

TYPED_TEST(ParticlesExportingTest, select_particles_into_flat_shifted)
{
    using Flat_t = ParticleArray<ParticleArrayOptions{
        TestFixture::dim, LayoutMode::AoS, StorageMode::VECTOR, TypeParam::alloc_mode}>;
    using Box_t = TestFixture::Box_t;

    auto [src, ref] = this->fill(false);
    // across the upper ghost layer, as a patch ghost selection
    Box_t const box{this->domain.upper - 1, this->ghost_box.upper};
    auto shift = ConstArray<int, TestFixture::dim>(0);
    shift[0]   = -TestFixture::cells;

    Flat_t dst;
    select_particles<ParticleType::PatchGhost>(src, dst, box, Point{shift});

    auto expected = this->keep_if(ref, [&](auto const& p) { return isIn(p, box); });
    for (auto& p : expected)
        p.iCell()[0] += shift[0];
    this->expect_same(dst, expected, "select shifted");
    EXPECT_EQ(count_particles<ParticleType::PatchGhost>(src, box), dst.size());
}

TYPED_TEST(ParticlesExportingTest, select_particles_into_same_layout)
{
    using Box_t = TestFixture::Box_t;

    auto [src, ref] = this->fill(false);
    Box_t const box{this->domain.lower, this->domain.lower + 4};
    auto dst = this->make_array();
    select_particles<ParticleType::Domain>(src, dst, box);
    if constexpr (requires { dst.template on_appended<ParticleType::Domain>(); })
        dst.template on_appended<ParticleType::Domain>();
    this->expect_same(dst, this->keep_if(ref, [&](auto const& p) { return isIn(p, box); }),
                      "select");
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
