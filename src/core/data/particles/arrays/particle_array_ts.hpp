#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP


#include "core/def.hpp"
#include "core/def/phare_config.hpp" // IWYU pragma: keep

#include "core/operators.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"

#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_translation_tracking.hpp"

namespace PHARE::core
{

template<typename Particles>
class ParticlesTile : public Box<std::int32_t, Particles::dimension>
{
    auto constexpr static dim = Particles::dimension;
    using This                = ParticlesTile<Particles>;
    using Super               = Box<std::int32_t, dim>;

    template<typename P>
    friend class ParticlesTile;


public:
    ParticlesTile(Super const& box, std::size_t const ghost_cells)
        : Super{box}
        , particles{make_particles<Particles>(box, ghost_cells)}
    {
    }

    ParticlesTile(ParticlesTile const&)            = default;
    ParticlesTile(ParticlesTile&&)                 = default;
    ParticlesTile& operator=(ParticlesTile const&) = default;
    ParticlesTile& operator=(ParticlesTile&&)      = default;

    template<typename Ps>
    ParticlesTile(ParticlesTile<Ps>& tile)
        : Super{tile}
        , particles{tile()}
    {
    }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

    auto& operator()() { return particles; }
    auto& operator()() const { return particles; }

    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }


    template<typename P>
    void reset(ParticlesTile<P>& tile)
    {
        (*this)().reset(tile());
    }

private:
    Particles particles;

    std::array<ParticlesTile*, 7> _links = ConstArray<ParticlesTile*, 7>(nullptr);
};

template<typename Particles>
struct CrossTileCopyDAO;

template<typename Particles>
class TileSetSpan : public ParticleTranslationTrackerSpan<Particles::dimension>
{
    using Base = ParticleTranslationTrackerSpan<Particles::dimension>;

    template<typename>
    friend struct CrossTileCopyDAO;

protected:
    using SIZE_T = default_span_size_t;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;

    using This               = TileSetSpan<Particles>;
    using lobox_t            = Box<std::uint32_t, dim>;
    using per_tile_particles = Particles;

    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::size;

    using Base::add_into_;
    using Base::box_;
    using Base::cap_;
    using Base::gap_idx_;
    using Base::gaps_;
    using Base::ghost_box_;
    using Base::left_;
    using Base::size_;

private:
    using locell_t = std::array<std::uint32_t, dim>;

    template<typename TileSetArray>
    auto resolve(TileSetArray& arr)
    {
        if constexpr (TileSetArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
            return arr.particles_views_.make_view();
    }

    template<typename TileSetArray>
    auto resolve_gaps(TileSetArray& arr)
    {
        if constexpr (TileSetArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
            return *arr.gap_views_;
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;

    template<typename TileSetArray>
    TileSetSpan(TileSetArray& arr)
        : Base{{.gaps            = resolve_gaps(arr),
                .gap_idx         = arr.gap_idx_,
                .add_into        = arr.add_into_,
                .cap             = arr.cap_,
                .left            = arr.left_,
                .size            = arr.size(),
                .box             = arr.box_,
                .ghost_box       = arr.ghost_box(),
                .local_ghost_box = arr.local_box()}}
        , particles_{resolve(arr)}
        , cell_size_{arr.cell_size_}
    {
    }

    auto size(std::size_t const& idx) const { return particles_.data()[idx].size(); }
    void resize(std::size_t s) { particles_.s = s; }


    template<typename TileSetArray>
    void reset(TileSetArray& arr)
    {
        particles_.reset(arr.particles_);
        size_ = arr.size();
    }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return (*particles_.at(cell))(); }
    auto& operator()(locell_t const& cell) const { return (*particles_.at(cell))(); }


    template<auto type, typename... Args>
    void sync(Args&&... args);


    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        size_ = 0;
    }


    auto local_tile_cell(std::array<int, dim> const& cell) const
    {
        PHARE_ASSERT(particles_.at(local_cell(cell)));
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

protected:
    void sync_tile_add_new(std::size_t const tidx);
    void sync_tile_rm_left(std::size_t const tidx);


    void static sort(auto from, auto to) { std::sort(from, to); }


    TileSetView<ParticlesTile<Particles>> particles_;
    NdArrayView<dim, SIZE_T> cell_size_;

}; // TileSetSpan


template<typename Particles>
class TileSetVector : public ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>
{
    using This = TileSetVector<Particles>;
    using Base = ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>;
    template<typename P>
    friend class TileSetSpan;

protected:
    using SIZE_T = default_span_size_t;

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
    using PSpan_t                      = Particles::view_t;
    using per_tile_particles           = Particles;
    using SpnTile                      = ParticlesTile<PSpan_t>;
    using VecTile                      = ParticlesTile<Particles>;
    using size_t_vector                = std::vector<std::size_t>;

    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::zero_bookkeeping;

    using Base::add_into_;
    using Base::box_;
    using Base::cap_;
    using Base::cell_size_;
    using Base::gap_idx_;
    using Base::gap_views_;
    using Base::gaps_;
    using Base::ghost_box_;
    using Base::ghost_cells_;
    using Base::left_;
    using Base::total_size;

    TileSetVector(box_t const& box, auto const ghost_cells)
        : Base{box, ghost_cells}
    {
        reset_views();
        TileSet<VecTile, alloc_mode>::build_links(particles_);
        TileSet<SpnTile, alloc_mode>::build_links(particles_views_);
    }

    TileSetVector(TileSetVector&& that)
        : Base{that.box_, that.ghost_cells_}
        , particles_{std::move(that.particles_)}
    {
        on_appended(); // without std::swap does not work well
    }

    TileSetVector(TileSetVector const& that)
        : Base{that.box_, that.ghost_cells_}
        , particles_{that.particles_.copy(TileSetter<dim>{that.box_, that.ghost_cells_})}
    {
        on_appended();
    }

    // resuse constructors rather than duplicate
    TileSetVector& operator=(TileSetVector&& that)
    {
        if (this == &that)
            return *this;
        this->~TileSetVector();
        new (this) TileSetVector(std::move(that));
        return *this;
    }
    TileSetVector& operator=(TileSetVector const& that)
    {
        if (this == &that)
            return *this;
        this->~TileSetVector();
        new (this) TileSetVector(that);
        return *this;
    }

    auto size() const { return total_size; }
    auto size(std::array<std::uint32_t, dim> const& icell) const { return cell_size_(icell); }
    auto size(std::size_t const& idx) const { return cell_size_.data()[idx]; }

    template<typename Iterator>
    auto erase(Iterator first, Iterator last)
    {
        return particles_.erase(particles_.begin() + first.curr_pos,
                                particles_.begin() + last.curr_pos);
    }

    void _inc(locell_t const& locell)
    {
        ++total_size;
        ++cell_size_(locell);
    }

    template<bool inc_ = true>
    void emplace_back(Particle_t const& p)
    {
        auto const locell = local_cell(p.iCell());
        assert(particles_.at(locell));
        (*particles_.at(locell))().emplace_back(p);
        if constexpr (inc_)
            _inc(locell);
    }

    template<bool inc_ = true>
    void emplace_back(Particles& dst, Particles const& src, std::size_t const& idx)
    {
        dst.emplace_back(src, idx);
    }

    template<typename... Args>
    void emplace_back(double const weight, Args&&... args)
    {
        this->emplace_back(Particle_t{weight, args...});
    }

    template<typename V>
    static auto& get_vec(V& v)
    {
        return v;
    }

    void push_back(Particle_t&& p) { emplace_back(p); }
    void push_back(Particle_t const& p) { emplace_back(p); }


    void reset_views()
    {
        update_from(reset_particle_views_fn(), particles_views_);
        update_from(reset_gap_views_fn(), gap_views_);
    }


    auto reset_gap_views_fn()
    {
        return [&](auto const i) { return make_span(*(gaps_.data() + i)); };
    }
    auto reset_particle_views_fn()
    {
        return [&](std::size_t const i) { return SpnTile{particles_[i]}; };
    }



    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    auto& views() { return particles_views_; }
    auto& views() const { return particles_views_; }
    auto& view(locell_t const& cell) { return (*particles_views_.at(cell))(); }

    // refreshes recount/gap-sizing/views after particles were added without going through
    // move_check (a raw append) - see on_moved() for the post-move_check pipeline.
    template<auto type = ParticleType::Domain>
    void on_appended();

    // applies move_check-registered moves: realloc/resize, span-side cross-tile copy, then
    // on_appended() to finalize.
    template<auto type>
    void on_moved(auto&&... args);

    template<auto type>
    void sync_moved(); // realloc ahead of the span-side copy (see sync_check_realloc)

    template<auto type>
    void trim();


    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        for (auto& tile : particles_views_)
            tile().clear();

        // cell_size_ here is cached per-tile (set in on_appended()), not read back from
        // each tile's own size() - reset it and the rest of the move_check bookkeeping
        // the same way the constructor does, or size(icell) stays stale after clear().
        zero_bookkeeping();

        reset_views();
        total_size = 0;
    }


    void static resize(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        resize(ps.particles_, s, copy);
    }

    template<typename V>
    void static resize(V& v, std::size_t const& s, bool const& copy = true)
    {
        v.resize(s);
    }

    void static reserve(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        reserve(ps.particles_, s, copy);
    }

    template<typename V>
    void static reserve(V& v, std::size_t const& s, bool const& copy = true)
    {
        v.reserve(s);
    }




    auto local_tile_cell(std::array<int, dim> const& cell) const
    {
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }



protected:
    template<auto type>
    void sync_check_realloc();

    auto on_tiles(auto&& fn)
    {
        for (auto& tile : particles_)
            fn(tile);
    }


    // tiles build from amrbox, but `.at()` function maps ghosts box
    TileSet<VecTile, alloc_mode> particles_{TileSetter<dim>{box_, ghost_cells_}, ghost_cells_};
    TileSet<SpnTile, alloc_mode> particles_views_ = TileSet<SpnTile, alloc_mode>::make_from(
        [](auto& tile) -> auto& { return tile; }, TileSetter<dim>{box_, ghost_cells_}, particles_);

}; // TileSetVector<Particles>




template<typename Particles>
template<auto type>
void TileSetVector<Particles>::trim() // change to erase(box)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
}



template<typename Particles>
template<auto type>
void TileSetVector<Particles>::sync_check_realloc()
{
    PHARE_LOG_SCOPE(3, "TileSetVector::sync_check_realloc");

    using enum LayoutMode;

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        auto const lix       = local_cell(particles_[i].lower);
        auto const& nu       = add_into_(lix);
        auto& real           = particles_[i]();
        auto const& left     = gap_idx_(lix);
        auto const& old_size = real.size();
        cell_size_(lix)      = old_size;
        auto const& new_size = real.size() + nu - left;
        reserve(real, real.size() + nu); // we must add new before removing leavers
        resize(real, new_size);
        cap_(lix) = real.capacity();
        assert(real.size() == cell_size_(lix) + nu - left);
        assert(cap_(lix) >= real.size());
        add_into_(lix) = 0;
    };
}


template<typename Particles>
template<auto type>
void TileSetVector<Particles>::sync_moved()
{
    sync_check_realloc<type>();
    reset_views();

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        auto const lix = local_cell(particles_[i].lower);
        particles_views_[i]().resize(cell_size_(lix));
    };
}

template<typename Particles>
template<auto type>
void TileSetVector<Particles>::on_moved(auto&&... args)
{
    sync_moved<type>();                                       // realloc + resize
    TileSetSpan<PSpan_t>{*this}.template sync<type>(args...); // cross-tile copy

    on_appended<type>(); // finalize: recount, size gaps, reset views
}

template<typename Particles>
template<auto type>
void TileSetVector<Particles>::on_appended()
{
    static_assert(all_are<ParticleType>(type));
    static_assert(type != ParticleType::All);
    PHARE_LOG_SCOPE(3, "TileSetVector::on_appended");

    total_size = 0;
    for (auto const& tile : particles_)
    {
        total_size += tile().size();
        cell_size_(local_cell(tile.lower)) = tile().size();
    }

    using enum LayoutMode;

    on_tiles([&](auto& tile) {
        auto const lix  = local_cell(tile.lower);
        auto const& cs  = cell_size_(lix);
        auto const& cap = tile().capacity();
        auto& gaps      = gaps_(lix);
        if (gaps.size() < cs)
        {
            reserve(gaps, cap, false);
            resize(gaps, cs, false);
        }
        cap_(lix) = cap;
        if constexpr (any_in(layout_mode, LayoutMode::AoSMapped))
            tile().remap();
    });

    PHARE_LOG_SCOPE(3, "TileSetVector::sync::reset_views");
    reset_views();
}


template<typename Super_>
struct TileSetParticles : public Super_
{
    using Super              = Super_;
    using This               = TileSetParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_tile_particles = typename Super::per_tile_particles;

    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::local_box;
    using Super::particles_;
    using Super::size;


    TileSetParticles(TileSetParticles&& that)
        : Super{std::forward<Super>(that)}
    {
    }

    TileSetParticles& operator=(TileSetParticles&&)      = default;
    TileSetParticles(TileSetParticles const&)            = default;
    TileSetParticles& operator=(TileSetParticles const&) = default;

    template<typename... Args>
    TileSetParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    // no begin()/end() here - AoSTS is tiled, so top-level particle iteration is not
    // supported (see ParticleArray::begin()/end() in particle_array.hpp); use
    // enumerate_tiles()/enumerate()/per_particle() instead.

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }


    template<auto particle_type>
    auto& move_check(auto const& pt, std::size_t const idx, auto& particle)
    {
        static_assert(particle_type == ParticleType::Domain);
        auto const& newcell = particle.iCell();
        if (array_equals(Super::local_tile_cell(newcell), pt.tile_cell))
            return *this;

        bool constexpr static ATOMIC = true;
        using Op                     = Operators<typename Super::SIZE_T, ATOMIC>;

        auto& gidx      = Super::gap_idx_(pt.tile_cell);
        auto const nidx = Op{gidx}.increment_return_old();
        auto& gaps      = Super::gaps_(pt.tile_cell);
        assert(nidx < gaps.size());
        gaps[nidx] = idx;

        using enum LayoutMode;
        if (isIn(newcell, Super::ghost_box()))
        {
            auto const nc = Super::local_tile_cell(newcell);
            Op{Super::add_into_(nc)}.increment_return_old();
        }

        return *this;
    }

    auto& evictor() {}


    template<typename T>
    struct index_wrapper;
    auto operator[](std::size_t const& s) { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const { return index_wrapper<This const>{this, s}; }

    void print() const {}
    void check() const;


    auto max_size() const
    {
        return max_from(this->particles_,
                        [](auto const& v, auto const& i) { return v.data()[i].size(); });
    }

    auto nbr_particles_in(std::array<int, dimension> const arr) const
    {
        auto const& tile = *particles_.at(Super::local_cell(arr));
        return tile().size() / tile.size(); // confusing? :)
    }

    auto nbr_particles_in(Box<int, dimension> const box) const
    {
        std::size_t n_particles = 0;

        traverse_tiles(this->particles_, box, [&](auto& tile) {
            auto overlap = box * grow(tile, this->ghost_cells_);
            assert(overlap);
            for (auto const bix : *overlap)
                n_particles += this->cell_size_(this->local_cell(bix));
        });

        return n_particles;
    }


}; // TileSetParticles<Super>


template<typename Particles>
template<auto type, typename... Args>
void TileSetSpan<Particles>::sync(Args&&... args)
{
    PHARE_LOG_SCOPE(3, "TileSetSpan::sync(stream)");

    for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
        sync_tile_add_new(tidx);

    for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
        sync_tile_rm_left(tidx);
}




template<typename Particles>
struct CrossTileCopyDAO
{
    bool constexpr static ATOMIC = true;
    using Op                     = Operators<default_span_size_t, ATOMIC>;
    using Tile                   = ParticlesTile<typename Particles::per_tile_particles>;

    Particles& ps;
    std::size_t src_tile_idx;
    Tile& tile = ps()[src_tile_idx];

    void copy_in()
    {
        auto& real         = tile();
        auto const bix     = ps.local_cell(tile.lower);
        auto const& n_gaps = ps.gap_idx_(bix);
        {
            auto& gaps = ps.gaps_(bix);
            ps.sort(gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
        }
        auto& left       = ps.left_(bix);
        auto const& gaps = ps.gaps_(bix);
        for (std::size_t i = 0; i < n_gaps; ++i)
        {
            auto const& gidx    = gaps[n_gaps - (1 + i)];
            auto const& newcell = ps.local_tile_cell(real.iCell(gidx));
            auto& ntile         = (*ps().at(newcell));
            auto& nparts        = ntile();
            auto const npidx    = next_index(nparts.size(), nparts.size_address());
            PHARE_ASSERT(npidx < ps.cap_(newcell));
            PHARE_ASSERT(not isIn(real[gidx], tile));
            PHARE_ASSERT(isIn(real[gidx], ntile));
            nparts.assign(real, gidx, npidx);
            PHARE_ASSERT(isIn(nparts[npidx], ntile));
            ++left;
        }
        PHARE_ASSERT(left == n_gaps);
    }

    void rm_left()
    {
        auto& real       = tile();
        auto const bix   = ps.local_cell(tile.lower);
        auto const& gaps = ps.gaps_(bix);
        auto& left       = ps.left_(bix);
        auto& gaps_size  = ps.gap_idx_(bix);
        while (left)
        {
            auto const& rsize = real.size();
            auto const& pidx  = gaps[gaps_size - 1];
            if (pidx != rsize - 1)
                real.assign(rsize - 1, pidx);
            real.pop_back();
            --gaps_size;
            --left;
        }
        PHARE_ASSERT(left == 0);
        PHARE_ASSERT(gaps_size == 0);
    }

    auto static next_index(auto npidx, auto const size_address)
    {
        while (true)
        {
            auto inc = npidx + 1;
            auto old = Op::compare_and_swap(size_address, npidx, inc);
            if (npidx != old)
            {
                ++npidx;
                continue;
            }
            else
                break;
        }
        return npidx;
    }
};



template<typename Particles>
void TileSetSpan<Particles>::sync_tile_add_new(std::size_t const tidx)
{
    CrossTileCopyDAO<std::decay_t<decltype(*this)>>{*this, tidx}.copy_in();
}


template<typename Particles>
void TileSetSpan<Particles>::sync_tile_rm_left(std::size_t const tidx)
{
    CrossTileCopyDAO<std::decay_t<decltype(*this)>>{*this, tidx}.rm_left();
}



template<typename Super>
void TileSetParticles<Super>::check() const
{
    // static constexpr auto alloc_mode  = Super::alloc_mode;
    using enum LayoutMode;
    if constexpr (storage_mode == StorageMode::VECTOR)
    {
        for (auto const& tile : particles_())
        {
            // auto const& tile = *particles_.at(gix);
            // assert(tile().capacity() > tile.size());

            auto const bix   = this->local_cell(tile.lower);
            auto const& gaps = this->gaps_(bix);
            auto& left       = this->left_(bix);
            auto& gaps_size  = this->gap_idx_(bix);

            assert(gaps_size == 0);
            assert(left == 0);
            assert(gaps.size() >= tile().size());
        }
    }
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP */
