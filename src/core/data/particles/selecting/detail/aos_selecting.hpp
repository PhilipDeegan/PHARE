
#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING


// #include "core/utilities/memory.hpp"
#include "core/logger.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"


#include <stdexcept>

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


//

// AoS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    for (auto const& p : src)
        if (isIn(p, box))
            dst.emplace_back(p);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    for (auto const& p : src)
        if (isIn(p, box))
            dst.emplace_back(shift_particle(p, shifter));
}

//

// AoSTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    // Domain particles live in a tile's own (domain) box; patch/level ghost particles
    // live only in a tile's own (clamped, uniquely-owned) ghost cells -- see
    // TileSetParticles::nbr_particles_in for the same grow()-based reasoning.
    auto const tile_box = [&](auto& tile) {
        if constexpr (particle_type == ParticleType::Domain)
            return *tile;
        else
        {
            auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
            return grow(tile, ghost_cells);
        }
    };

    for (auto& src_tile : const_cast<SrcParticles&>(src)())
        if (auto const src_box_overlap = tile_box(src_tile) * box)
            for (auto const& p : src_tile())
                if (isIn(p, *src_box_overlap))
                    dst.push_back(p);

    if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
        dst.template on_appended<particle_type>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    using enum LayoutMode;
    static_assert(any_in(DstParticles::layout_mode, AoS, AoSTS), "Otherwise not implemented!");

    auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        if (auto const src_overlap = grow(src_tile, ghost_cells) * box)
        {
            if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
            {
                for (auto const& p : src_tile())
                    if (isIn(p, *src_overlap))
                        dst.push_back(shift_particle(p, shifter));
            }
            else
            {
                auto const copy = [&](auto& s, auto& d, auto const& b) {
                    auto const in_dst = partition_particles(s, b);
                    if (in_dst.size() == 0)
                        return;

                    auto const old_dst_size = d.size();
                    d.resize(d.size() + in_dst.size());

                    mem::copy<CPU>(d.data() + old_dst_size, s.data(), in_dst.size());
                    if (any(shifter, [](auto const e) { return e != 0; }))
                        for (std::size_t i = old_dst_size; i < d.size(); i++)
                            d[i].iCell() = *((shifter) + d[i].iCell());
                };
                copy(src_tile(), dst, *src_overlap);
            }
        }
    }

    if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
        dst.template on_appended<particle_type>();
}


//

// AoSPCTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPCTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    // ghost_box(), not the tile's bare domain box, for ghost-layer selections: a
    // tile's per-cell buckets extend into its own (clamped, uniquely-owned) ghost
    // cells, which is where patchGhost/levelGhost particles actually live -- see
    // TileSetParticles::nbr_particles_in
    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        auto& src_pc = src_tile();
        auto const tile_box
            = particle_type == ParticleType::Domain ? src_pc.box() : src_pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_src_box = src_pc.local_box(*overlap);
            for_each_point(lcl_src_box, [&](auto const& cell) {
                for (auto const& p : src_pc(cell))
                    dst.push_back(p);
            });
        }
    }

    // finalize any destination with per-cell/per-tile bookkeeping (AoSPC, AoSTS, AoSPCTS)
    if constexpr (requires { dst.template on_appended<particle_type>(); })
        dst.template on_appended<particle_type>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPCTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        auto& src_pc = src_tile();
        auto const tile_box
            = particle_type == ParticleType::Domain ? src_pc.box() : src_pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_src_box = src_pc.local_box(*overlap);
            for_each_point(lcl_src_box, [&](auto const& cell) {
                for (auto const& p : src_pc(cell))
                    dst.push_back(shift_particle(p, shifter));
            });
        }
    }

    // finalize any destination with per-cell/per-tile bookkeeping (AoSPC, AoSTS, AoSPCTS)
    if constexpr (requires { dst.template on_appended<particle_type>(); })
        dst.template on_appended<particle_type>();
}


//

// AoSCMTS, CPU
// Domain particles live only in their own tile, ghost-layer particles (level ghosts) only in
// their clamp-owner tile (TileSet::tag_cells_), so ghost cells are taken from that tile.
// fn(tile, box) gets each tile's owned sub-box of box, to use with the tile's cellmap
template<ParticleType particle_type, typename SrcParticles, typename Fn>
void on_owned_cell_mapped_boxes(SrcParticles const& src, auto const& box, Fn&& fn)
{
    using Box_t       = Box<int, SrcParticles::dimension>;
    auto const& tiles = src();
    for (auto const& tile : tiles)
    {
        if constexpr (particle_type == ParticleType::Domain)
        {
            if (auto const overlap = *tile * box)
                fn(tile, *overlap);
        }
        else if (auto const overlap = grow(*tile, src.ghost_cells()) * box)
            src.amr_lcl_idx(*overlap, [&](auto const& bix, auto const& lix) {
                if (tiles.at(lix) == &tile)
                    fn(tile, Box_t{bix, bix});
            });
    }
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSCMTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    on_owned_cell_mapped_boxes<particle_type>(
        src, box, [&](auto const& tile, auto const& b) { tile().export_particles(b, dst); });
    dst.template on_appended<particle_type>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSCMTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    auto const offseter = [&](auto const& particle) { return shift_particle(particle, shifter); };
    on_owned_cell_mapped_boxes<particle_type>(src, box, [&](auto const& tile, auto const& b) {
        tile().export_particles(b, dst, offseter);
    });
    dst.template on_appended<particle_type>();
}


//

// AoS, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

// AoSPC, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    for_each_point_zip(lcl_src_box, lcl_dst_box, [&](auto const& src_cell, auto const& dst_cell) {
        auto& sv = src(src_cell);
        auto& dv = dst(dst_cell);
        dv.reserve(dv.size() + sv.size());
        for (auto const& p : sv)
            dv.emplace_back(p);
    });
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - fn); // BEWARE
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    for_each_point_zip(lcl_src_box, lcl_dst_box, [&](auto const& src_cell, auto const& dst_cell) {
        auto& sv = src(src_cell);
        auto& dv = dst(dst_cell);
        dv.reserve(dv.size() + sv.size());
        for (auto p : sv)
        {
            p.iCell() = (Point{p.iCell()} + fn).toArray();
            dv.emplace_back(p);
        }
    });
}

//

// AoSPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    for_each_point_zip(lcl_src_box, lcl_dst_box, [&](auto const& src_cell, auto const& dst_cell) {
        auto& sv = src(src_cell);
        auto& dv = dst(dst_cell);
        dv.reserve(dv.size() + sv.size());
        for (auto const& p : sv)
            dv.emplace_back(p);
    });
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shift)
{ // unoptimized
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - shift);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    for_each_point_zip(lcl_src_box, lcl_dst_box, [&](auto const& src_cell, auto const& dst_cell) {
        auto& sv = src(src_cell);
        auto& dv = dst(dst_cell);
        dv.reserve(dv.size() + sv.size());
        for (auto const& pit : sv)
        {
            auto p    = pit.copy();
            p.iCell() = (Point{p.iCell()} + shift).toArray();
            dv.emplace_back(p);
        }
    });
}

//

// AoSMapped, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    src.export_particles(box, dst);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    auto const offseter = [&](auto const& particle) { return shift_particle(particle, shifter); };

    src.export_particles(box, dst, offseter);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSMapped, CPU>::count(SrcParticles const& src, box_t const& box)
{
    return src.nbr_particles_in(box);
}

//

// AoS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto const& p : src)
        if (isIn(p, box))
            ++n;
    return n;
}

// AoS, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoS, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    throw std::runtime_error("finish this");
    return 0;
}

// AoSTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto& tile : const_cast<SrcParticles&>(src)())
    {
        auto const box_for_overlap = [&]() {
            if constexpr (particle_type == ParticleType::Domain)
                return *tile;
            else
            {
                auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
                return grow(tile, ghost_cells);
            }
        }();
        if (auto const overlap = box_for_overlap * box)
            for (auto const& p : tile())
                if (isIn(p, *overlap))
                    ++n;
    }
    return n;
}

// AoSPCTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPCTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto& tile : const_cast<SrcParticles&>(src)())
    {
        auto& pc            = tile();
        auto const tile_box = particle_type == ParticleType::Domain ? pc.box() : pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_box = pc.local_box(*overlap);
            for_each_point(lcl_box, [&](auto const& cell) { n += pc(cell).size(); });
        }
    }
    return n;
}

// AoSCMTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSCMTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    on_owned_cell_mapped_boxes<particle_type>(
        src, box, [&](auto const& tile, auto const& b) { n += tile().nbr_particles_in(b); });
    return n;
}

// AoSPC, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPC, CPU>::count(SrcParticles const& src, box_t const& box)
{
    auto const overlap = src.ghost_box() * box; // box may extend past src's own range
    if (not overlap)
        return 0;
    auto const lcl_box = src.local_box(*overlap);
    std::size_t n      = 0;
    for_each_point(lcl_box, [&](auto const& cell) { n += src(cell).size(); });
    return n;
}

// AoSPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPC, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    auto const lcl_box = src.local_box(box);
    std::size_t n      = 0;
    for_each_point(lcl_box, [&](auto const& cell) { n += src(cell).size(); });
    return n;
}


//

// AoSTS / AoSPCTS, GPU_UNIFIED
// the CPU selectors' tile / box overlaps, then on the device: one block per overlapping tile
// counts the selected particles, then gathers them (shifted) into managed memory, which is
// handed to dst: appended on the device for AoSTS GPU, copied in bulk for flat arrays (e.g. the
// pack buffer). Only particle positions decide, so the order differs from the CPU selectors

namespace detail::gpu_select
{
    template<std::size_t dim>
    struct TileSelection
    {
        std::size_t tile;
        Box<int, dim> box;               // AMR, the selected part of the tile's particles
        Box<std::uint32_t, dim> lcl_box; // AoSPCTS: box in the tile's per cell indexing
    };

    // as the CPU selectors: Domain particles live in a tile's own box, ghost particles in its
    //  (clamped, uniquely owned) ghost cells
    template<ParticleType type, typename Src, typename Box_t>
    auto tile_selections(Src const& src, Box_t const& box)
    {
        auto constexpr dim = Src::dimension;
        std::vector<TileSelection<dim>> sels;
        auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
        for (std::size_t ti = 0; ti < src().size(); ++ti)
        {
            auto const& tile = src()[ti];
            auto const tbox  = type == ParticleType::Domain ? Box<int, dim>{*tile}
                                                            : grow(*tile, ghost_cells);
            if (auto const overlap = tbox * box)
            {
                if constexpr (Src::layout_mode == LayoutMode::AoSPCTS)
                    sels.push_back({ti, *overlap, tile().local_box(*overlap)});
                else
                    sels.push_back({ti, *overlap, {}});
            }
        }
        return sels;
    }

    // block threads: fn(particle) for every selected particle of the selection's tile
    template<typename Launcher, typename View>
    void on_selected(View& view, auto const& sel, auto&& fn) _PHARE_ALL_FN_
    {
        if constexpr (View::layout_mode == LayoutMode::AoSPCTS)
        {
            auto& cps = view()[sel.tile]();
            for (std::size_t c = Launcher::thread_idx(); c < sel.lcl_box.size();
                 c += Launcher::threads())
            {
                auto& parts = cps(point_at(sel.lcl_box, c));
                for (std::size_t i = 0; i < parts.size(); ++i)
                    fn(parts[i]);
            }
        }
        else
        {
            auto& parts = view()[sel.tile]();
            for (std::size_t i = Launcher::thread_idx(); i < parts.size(); i += Launcher::threads())
                if (isIn(parts[i], sel.box))
                    fn(parts[i]);
        }
    }

#if PHARE_HAVE_MKN_GPU
    // the selected particles of src in box, shifted, gathered into managed memory
    template<ParticleType type, typename Src, typename Box_t>
    auto gather(Src const& src, Box_t const& box, auto const& shift)
    {
        auto constexpr dim = Src::dimension;
        using Particle_t   = typename Src::Particle_t;
        using Launcher     = gpu::BlockPerLauncher<true>;
        using Op           = Operators<unsigned long long, true, true>;

        gpu::Vec_t<Particle_t> out;
        auto const host_sels = tile_selections<type>(src, box);
        if (host_sels.empty())
            return out;

        gpu::Vec_t<TileSelection<dim>> const sels(host_sels.begin(), host_sels.end());
        gpu::Vec_t<unsigned long long> counter(1, 0);
        auto const* const sp = sels.data();
        auto* const cnt      = counter.data();
        auto view            = *const_cast<Src&>(src);

        Launcher{sels.size()}([=] _PHARE_ALL_FN_() mutable {
            on_selected<Launcher>(view, sp[Launcher::block_idx()],
                                  [&](auto const&) { Op{*cnt}.increment_return_old(); });
        });

        out.resize(counter[0]);
        if (out.empty())
            return out;
        counter[0]          = 0;
        auto* const o       = out.data();
        auto const shift_by = for_N_make_array<dim>([&](auto i) { return int(shift[i]); });

        Launcher{sels.size()}([=] _PHARE_ALL_FN_() mutable {
            on_selected<Launcher>(view, sp[Launcher::block_idx()], [&](auto const& p) {
                auto& q = o[Op{*cnt}.increment_return_old()];
                q       = p;
                for (std::size_t d = 0; d < dim; ++d)
                    q.iCell()[d] += shift_by[d];
            });
        });
        return out;
    }
#endif // PHARE_HAVE_MKN_GPU

    template<ParticleType type, typename Dst>
    void hand_to(Dst& dst, auto const& staged)
    {
        if constexpr (Dst::layout_mode == LayoutMode::AoSTS
                      and Dst::alloc_mode == AllocatorMode::GPU_UNIFIED)
            dst.template append_on_device<type>(staged.data(), staged.size());
        else
        {
            if constexpr (Dst::layout_mode == LayoutMode::AoS)
            {
                auto const old = dst.size();
                dst.resize(old + staged.size());
                std::copy(staged.begin(), staged.end(), dst.data() + old);
            }
            else
                for (auto const& p : staged)
                    dst.push_back(p);
            if constexpr (requires { dst.template on_appended<type>(); })
                dst.template on_appended<type>();
        }
    }

    template<ParticleType type, typename Src, typename Dst, typename Box_t>
    void select(Src const& src, Dst& dst, Box_t const& box, auto const& shift)
    {
#if PHARE_HAVE_MKN_GPU
        hand_to<type>(dst, gather<type>(src, box, shift));
#else
        throw std::runtime_error("GPU particle selection needs mkn.gpu");
#endif
    }

    template<typename Src>
    auto no_shift()
    {
        return ConstArray<int, Src::dimension>(0);
    }

} // namespace detail::gpu_select


template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSTS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    using namespace detail::gpu_select;
    detail::gpu_select::select<particle_type>(src, dst, box, no_shift<SrcParticles>());
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSTS, GPU_UNIFIED>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    detail::gpu_select::select<particle_type>(src, dst, box, shifter);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSTS, GPU_UNIFIED>::count(SrcParticles const& src,
                                                         box_t const& box)
{
#if PHARE_HAVE_MKN_GPU
    using namespace detail::gpu_select;
    using Launcher = gpu::BlockPerLauncher<true>;
    using Op       = Operators<unsigned long long, true, true>;

    auto const host_sels = tile_selections<particle_type>(src, box);
    if (host_sels.empty())
        return 0;
    gpu::Vec_t<TileSelection<SrcParticles::dimension>> const sels(host_sels.begin(),
                                                                  host_sels.end());
    gpu::Vec_t<unsigned long long> counter(1, 0);
    auto const* const sp = sels.data();
    auto* const cnt      = counter.data();
    auto view            = *const_cast<SrcParticles&>(src);
    Launcher{sels.size()}([=] _PHARE_ALL_FN_() mutable {
        on_selected<Launcher>(view, sp[Launcher::block_idx()],
                              [&](auto const&) { Op{*cnt}.increment_return_old(); });
    });
    return counter[0];
#else
    throw std::runtime_error("GPU particle counting needs mkn.gpu");
#endif
}


template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPCTS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    using namespace detail::gpu_select;
    detail::gpu_select::select<particle_type>(src, dst, box, no_shift<SrcParticles>());
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPCTS, GPU_UNIFIED>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    detail::gpu_select::select<particle_type>(src, dst, box, shifter);
}

// per cell sizes only, no particle data is read: the CPU version
template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPCTS, GPU_UNIFIED>::count(SrcParticles const& src,
                                                           box_t const& box)
{
    return ParticlesSelector<AoSPCTS, CPU>::template count<particle_type>(src, box);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING */
