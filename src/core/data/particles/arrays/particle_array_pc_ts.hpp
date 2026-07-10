#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP


#include "core/def.hpp"
#include "core/operators.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_pc.hpp"



namespace PHARE::core
{


template<typename PerCellParticles>
class PCParticlesTile : public Box<std::int32_t, PerCellParticles::dimension>
{
    static constexpr auto dim = PerCellParticles::dimension;
    using Super               = Box<std::int32_t, dim>;
    using This                = PCParticlesTile<PerCellParticles>;

public:
    static constexpr auto dimension = dim;

    PCParticlesTile(Super const& box, std::size_t const ghost_cells)
        : Super{box}
        , particles{box, ghost_cells}
    {
    }

    PCParticlesTile(PCParticlesTile const&)            = default;
    PCParticlesTile(PCParticlesTile&&)                 = default;
    PCParticlesTile& operator=(PCParticlesTile const&) = default;
    PCParticlesTile& operator=(PCParticlesTile&&)      = default;

    template<typename OtherPerCell>
    PCParticlesTile(PCParticlesTile<OtherPerCell>& other)
        : Super{other}
        , particles{other.particles}
    {
    }

    auto& operator()() _PHARE_ALL_FN_ { return particles; }
    auto& operator()() const _PHARE_ALL_FN_ { return particles; }
    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }

    PerCellParticles particles;

private:
    std::array<This*, 7> _links = ConstArray<This*, 7>(nullptr);
};


template<typename Particles>
struct PCCrossTileCopyDAO;


template<typename Particles>
class PCTileSetSpan
{
    friend struct PCTileSetParticlesService;
    template<typename>
    friend struct PCCrossTileCopyDAO;

protected:
    using SIZE_T = default_span_size_t;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;

    using This               = PCTileSetSpan<Particles>;
    using lobox_t            = Box<std::uint32_t, dim>;
    using per_tile_particles = Particles;
    using Tile_t             = PCParticlesTile<Particles>;

private:
    using locell_t = std::array<std::uint32_t, dim>;

    template<typename PCTileSetArray>
    auto resolve(PCTileSetArray& arr)
    {
        if constexpr (PCTileSetArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
            return arr.particles_views_.make_view();
    }

    template<typename PCTileSetArray>
    auto resolve_gaps(PCTileSetArray& arr)
    {
        if constexpr (PCTileSetArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
            return *arr.gap_views_;
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;

    template<typename PCTileSetArray>
    PCTileSetSpan(PCTileSetArray& arr)
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cell_size_{arr.cell_size_}
        , cap_{arr.cap_}
        , left_{arr.left_}
        , size_{arr.size()}
        , box_{arr.box_}
        , ghost_box_{arr.ghost_box_}
        , local_ghost_box_{arr.local_box()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    auto size(std::size_t const& idx) const _PHARE_ALL_FN_
    {
        return particles_.data()[idx]().size();
    }

    auto& box() const _PHARE_ALL_FN_ { return box_; }
    auto& ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }

    auto local_cell(std::array<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return local_cell(icell.toArray());
    }

    auto& local_box() const _PHARE_ALL_FN_ { return local_ghost_box_; }

    auto local_box(Box<int, dim> const& from) const _PHARE_ALL_FN_
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    auto& operator()() _PHARE_ALL_FN_ { return particles_; }
    auto& operator()() const _PHARE_ALL_FN_ { return particles_; }
    auto& operator()(locell_t const& cell) _PHARE_ALL_FN_ { return (*particles_.at(cell))(); }
    auto& operator()(locell_t const& cell) const _PHARE_ALL_FN_ { return (*particles_.at(cell))(); }

    auto local_tile_cell(std::array<int, dim> const& cell) const _PHARE_ALL_FN_
    {
        PHARE_ASSERT(particles_.at(local_cell(cell)));
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

    template<auto particle_type>
    auto& move_check(auto const& pt, std::size_t const idx, auto& particle) _PHARE_ALL_FN_
    {
        static_assert(particle_type == ParticleType::Domain);

        auto const& newcell = particle.iCell();

        if (array_equals(newcell, pt.icell))
            return *this; // old cell == new cell, no change required

        auto& old_tile           = *particles_.at(pt.tile_cell);
        auto& old_cell_particles = old_tile();

        if (isIn(newcell, box()) and not isIn(newcell, old_tile))
        {
            // tile change register — mirrors TileSetParticles::move_check: only register
            // the move here, the actual cross-tile copy happens later during sync.
            bool constexpr static ATOMIC = true;
            bool constexpr static GPU    = alloc_mode == AllocatorMode::GPU_UNIFIED;
            using Op                     = Operators<SIZE_T, ATOMIC, GPU>;

            auto const old_lcl_cell = local_cell(pt.icell);

            auto& gidx      = gap_idx_(old_lcl_cell);
            auto const nidx = Op{gidx}.increment_return_old();
            auto& gaps      = gaps_(old_lcl_cell);
            assert(nidx < gaps.size());
            gaps[nidx] = idx;

            Op{add_into_(local_cell(newcell))}.increment_return_old();
            return *this;
        }

        // cell change register: still owned by (or ghosted into) the current tile.
        // icell_changer's own ghost_box() check covers both "outside the domain, assumed
        // still in this tile's ghost box" and "inside the domain, inside this tile".
        old_cell_particles.icell_changer(particle, old_cell_particles.local_cell(pt.icell), idx,
                                         newcell);

        return *this;
    }

    template<auto type, typename... Args>
    void sync(Args&&... /*args*/) _PHARE_ALL_FN_
    {
    } // TODO

    void clear() {} // TODO

protected:
    void sync_tile_add_new(std::size_t const /*tidx*/) _PHARE_ALL_FN_ {} // TODO
    void sync_tile_rm_left(std::size_t const /*tidx*/) _PHARE_ALL_FN_ {} // TODO

    TileSetView<Tile_t> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> gap_idx_, add_into_, cell_size_, cap_, left_;
    std::size_t size_;
    Box<int, dim> box_, ghost_box_;
    lobox_t local_ghost_box_;

}; // PCTileSetSpan


template<typename Particles>
class PCTileSetVector
{
    using This                    = PCTileSetVector<Particles>;
    using SIZE_T                  = default_span_size_t;
    bool static constexpr c_order = true;

    template<typename P>
    friend class PCTileSetSpan;

    friend struct PCTileSetParticlesService;

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr layout_mode  = Particles::layout_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = ParticleDefaults<dim>::Particle_t;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::view_t;
    using per_tile_particles           = Particles;
    using Tile_t                       = PCParticlesTile<Particles>;
    using SpnTile                      = PCParticlesTile<PSpan_t>;

    template<typename T>
    using nd_array_t = NdArrayVector<dim, T, c_order, alloc_mode>;

    std::uint8_t constexpr static alloc_impl()
    {
        if (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
            return 1;
        return 1;
    }

    template<typename T>
    using vec_helper = PHARE::Vector<T, alloc_mode, alloc_impl()>;

    using size_t_vector = typename vec_helper<std::size_t>::vector_t;

    PCTileSetVector(box_t const& box, auto const ghost_cells)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
        , particles_{TileSetter<dim>{box, ghost_cells}, ghost_cells}
        , particles_views_{TileSet<SpnTile, alloc_mode>::make_from(
              [](auto& tile) -> auto& { return tile; }, TileSetter<dim>{box, ghost_cells},
              particles_)}
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        left_.zero();
        cap_.zero();

        TileSet<Tile_t, alloc_mode>::build_links(particles_);
        TileSet<SpnTile, alloc_mode>::build_links(particles_views_);
    }

    PCTileSetVector(PCTileSetVector&& that)            = default;
    PCTileSetVector(PCTileSetVector const& that)       = default;
    PCTileSetVector& operator=(PCTileSetVector&&)      = default;
    PCTileSetVector& operator=(PCTileSetVector const&) = default;

    auto size() const { return total_size; }
    auto size(locell_t const& icell) const { return cell_size_(icell); }
    auto size(std::size_t const& idx) const { return cell_size_.data()[idx]; }

    template<auto type = ParticleType::Domain>
    auto& reserve_ppc(std::size_t const& ppc)
    {
        for (auto& tile : particles_)
            tile().template reserve_ppc<type>(ppc);
        return *this;
    }

    void emplace_back(Particle_t const& p)
    {
        auto* tile = particles_.at(Point<int, dim>{p.iCell()});
        assert(tile);
        (*tile)().emplace_back(p);
        ++total_size;
    }

    void push_back(Particle_t const& p) { emplace_back(p); }

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }
    auto local_box() const
    {
        return box_from_zero_to_upper_minus_one(
            ghost_box_.shape().template toArray<std::uint32_t>());
    }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync()
    {
        if constexpr (PHASE == 1)
            sync_check_realloc<type>();
        total_size = 0;
        for (auto const& tile : particles_)
            total_size += tile().size();
    }

    template<auto type>
    void sync_moved()
    { /* TODO */
    }

    template<auto type>
    void trim()
    { /* TODO */
    }

    void reset_views()
    {
        update_from([&](std::size_t const i) { return SpnTile{particles_[i]}; }, particles_views_);
        update_from([&](std::size_t const i) { return make_span(*(gaps_.data() + i)); },
                    gap_views_);
    }

    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        total_size = 0;
    }

protected:
    template<auto type>
    void sync_check_realloc()
    {
        // cross-tile realloc after moves — TODO: implement gap-based realloc
        total_size = 0;
        for (auto const& tile : particles_)
            total_size += tile().size();
    }

    std::size_t ghost_cells_;
    Box<int, dim> box_, ghost_box_;

    nd_array_t<size_t_vector> gaps_{local_box().shape()};
    nd_array_t<Span<std::size_t>> gap_views_{local_box().shape()};

    nd_array_t<SIZE_T> gap_idx_{local_box().shape()};
    nd_array_t<SIZE_T> add_into_{local_box().shape()};
    nd_array_t<SIZE_T> left_{local_box().shape()};
    nd_array_t<SIZE_T> cap_{local_box().shape()};
    nd_array_t<SIZE_T> cell_size_{local_box().shape()};

    TileSet<Tile_t, alloc_mode> particles_;
    TileSet<SpnTile, alloc_mode> particles_views_;
    std::size_t total_size = 0;

}; // PCTileSetVector<Particles>


template<typename Super_>
struct PCTileSetParticles : public Super_
{
    using Super              = Super_;
    using This               = PCTileSetParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_tile_particles = typename Super::per_tile_particles;

    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::size;

    PCTileSetParticles(PCTileSetParticles&&)                 = default;
    PCTileSetParticles& operator=(PCTileSetParticles&&)      = default;
    PCTileSetParticles(PCTileSetParticles const&)            = default;
    PCTileSetParticles& operator=(PCTileSetParticles const&) = default;

    template<typename... Args>
    PCTileSetParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    template<typename T>
    struct iterator_impl;

    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args)
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T&>{*t, args...};
    }
    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, true); }
    auto end() _PHARE_ALL_FN_ { return it(this, true); }

    auto data() const _PHARE_ALL_FN_ { return static_cast<Particle_t const*>(nullptr); } // TODO
    auto data() _PHARE_ALL_FN_ { return static_cast<Particle_t*>(nullptr); }             // TODO

    template<auto particle_type>
    auto& move_check(auto const& /*pt*/, std::size_t const /*idx*/,
                     auto& /*particle*/) _PHARE_ALL_FN_
    {
        return *this; /* TODO */
    }

    auto nbr_particles_in(std::array<int, dimension> const /*arr*/) const
    {
        return std::size_t{0}; /* TODO */
    }

    auto nbr_particles_in(Box<int, dimension> const /*box*/) const
    {
        return std::size_t{0}; /* TODO */
    }

    auto max_size() const { return std::size_t{0}; /* TODO */ }

    void print() const {}
    void check() const {}

    template<typename T>
    struct index_wrapper;
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return index_wrapper<This const>{this, s};
    }

}; // PCTileSetParticles<Super>


template<typename Particles>
struct pc_ts_iterator_base
{
    pc_ts_iterator_base(Particles ps)
        : particles{ps}
    {
    }

    Particles particles;
    std::size_t l = 0, i = 0;
};

template<auto layout_mode, typename Particles>
struct pc_ts_iterator_storage : public pc_ts_iterator_base<Particles>
{
    using Super = pc_ts_iterator_base<Particles>;

    pc_ts_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    void set() { /* TODO */ }
};

template<typename T>
struct pc_ts_iterator_super
{
    using per_tile_particles = typename std::decay_t<T>::per_tile_particles;
    using value_type         = pc_ts_iterator_storage<per_tile_particles::layout_mode, T>;
};
template<typename T>
using pc_ts_iterator_super_v = typename pc_ts_iterator_super<T>::value_type;

template<typename OuterSuper>
template<typename T>
struct PCTileSetParticles<OuterSuper>::iterator_impl : public pc_ts_iterator_super_v<T>
{
    // auto static constexpr dimension = OuterSuper::dimension;
    // using Super                     = pc_ts_iterator_super_v<T>;
    // using outer_type                = std::decay_t<T>;
    // using difference_type           = std::size_t;
    // using iterator_category         = std::forward_iterator_tag;
    // using Particle_t                = typename OuterSuper::Particle_t;
    // using value_type                = Particle_t;
    // using pointer                   = Particle_t*;
    // using reference                 = Particle_t&;
    // using Super::l, Super::i, Super::particles;

    // // l is used as a flat particle index; i unused until proper tile traversal is implemented
    // iterator_impl(T& particles_, bool end = false)
    //     : Super{particles_}
    // {
    //     if (end)
    //         l = particles.size(); // sentinel: flat index == total count
    //     // else l = 0 (begin)
    // }
    // iterator_impl(iterator_impl&&)                 = default;
    // iterator_impl(iterator_impl const&)            = default;
    // iterator_impl& operator=(iterator_impl&&)      = default;
    // iterator_impl& operator=(iterator_impl const&) = default;

    // auto operator*() const { return Particle_t{}; } // TODO: return actual particle
    // auto& operator++() { ++l; return *this; }
    // auto operator++(int) { auto copy = *this; ++(*this); return copy; }

    // auto operator==(iterator_impl const& that) const { return l == that.l; }
    // auto operator!=(iterator_impl const& that) const { return !(*this == that); }
    // auto operator<(iterator_impl const& that) const { return l < that.l; }

    // auto copy() const _PHARE_ALL_FN_ { return Particle_t{}; } // TODO
};


template<typename Particles>
struct PCCrossTileCopyDAO
{
    using Tile = PCParticlesTile<typename Particles::per_tile_particles>;

    Particles& ps;
    std::size_t src_tile_idx;

    void copy_in() _PHARE_ALL_FN_ { /* TODO */ }
    void rm_left() _PHARE_ALL_FN_ { /* TODO */ }
};


template<auto type>
void sync_pc_ts(auto& particles, auto&&... args)
{
    particles.template sync_moved<type>();     // realloc
    (*particles).template sync<type>(args...); // cross-tile copy

    PHARE_DEBUG_DO({
        auto& tiles = particles();
        auto& views = particles.views();
        for (std::size_t i = 0; i < particles().size(); ++i)
        {
            assert(views[i]().size() == tiles[i]().size());
        }
    })

    particles.template sync<0, type>(); // finalize
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP */
