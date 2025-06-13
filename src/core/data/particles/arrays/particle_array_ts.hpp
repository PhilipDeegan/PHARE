#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP


#include "core/def.hpp"
#include "core/def/detail/thrust.hpp"
#include "core/def/phare_config.hpp"
#include "core/operators.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"

#include "core/utilities/box/box.hpp"
#include "core/utilities/span.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"



#include <tuple>
#include <iterator>
#include <stdexcept>


namespace PHARE::core
{



template<typename Particles, typename NdArray_t>
class ParticlesTile : public Box<std::int32_t, Particles::dimension>
{
    auto constexpr static dim = Particles::dimension;
    using This                = ParticlesTile<Particles, NdArray_t>;
    using Super               = Box<std::int32_t, dim>;

    template<typename P, typename N>
    friend class ParticlesTile;

    auto constexpr static ndarray_builder = [](auto const& box) {
        auto const cells = *(box.shape().as_unsigned());
        if constexpr (Particles::storage_mode == StorageMode::VECTOR)
            return NdArray_t{cells}; // allocating
        else
            return NdArray_t{0, cells}; // span with nullptr, updated later
    };



public:
    ParticlesTile(ParticlesTile const&)            = default;
    ParticlesTile(ParticlesTile&&)                 = default;
    ParticlesTile& operator=(ParticlesTile const&) = default;
    ParticlesTile& operator=(ParticlesTile&&)      = default;


    ParticlesTile(Super const& box, std::size_t const ghost_cells)
        : Super{box}
        , particles{make_particles<Particles>(box, ghost_cells)}

    {
    }

    template<typename Ps, typename Nd, auto S0 = Particles::storage_mode,
             auto S1  = Ps::storage_mode,
             typename = std::enable_if_t<S0 == StorageMode::SPAN and S0 != S1>>
    ParticlesTile(ParticlesTile<Ps, Nd>& tile)
        : Super{tile}
        , particles{tile(), tile().size()}

    {
    }




    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

    auto& operator()() _PHARE_ALL_FN_ { return particles; }
    auto& operator()() const _PHARE_ALL_FN_ { return particles; }

    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }


    template<typename P, typename N>
    void reset(ParticlesTile<P, N>& tile)
    {
        (*this)().reset(tile());
    }

private:
    Particles particles;

    std::array<ParticlesTile*, 7> _links = ConstArray<ParticlesTile*, 7>(nullptr);
};



template<typename Particles, std::uint8_t impl_ = 0>
class TileSetSpan
{
protected:
    using SIZE_T = unsigned long long int;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;
    auto static constexpr impl_v     = impl_;
    using This                       = TileSetSpan<Particles, impl_v>;
    using lobox_t                    = Box<std::uint32_t, dim>;
    using per_tile_particles         = Particles;

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
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cap_{arr.cap_}
        , left_{arr.left_}
        , size_{arr.size()}
        , box_{arr.box_}
        , ghost_box_{arr.ghost_box()}
        , local_ghost_box_{arr.local_box()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    auto size(std::size_t const& idx) const _PHARE_ALL_FN_ { return particles_.data()[idx].size(); }

    void resize(std::size_t s) { particles_.s = s; }


    template<typename TileSetArray>
    void reset(TileSetArray& arr)
    {
        particles_.reset(arr.particles_);
        size_ = arr.size();
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


    template<typename... Args>
    void sync(Args&&... args) _PHARE_ALL_FN_;


    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        size_ = 0;
    }



    auto local_tile_cell(std::array<int, dim> const& cell) const _PHARE_ALL_FN_
    {
        PHARE_ASSERT(particles_.at(local_cell(cell)));

        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

protected:
    void sync_tile_add_new(std::size_t const tidx) _PHARE_ALL_FN_;
    void sync_tile_rm_left(std::size_t const tidx) _PHARE_ALL_FN_;

    void static sort(auto from, auto to) _PHARE_ALL_FN_
    {
        int sorted = 0;
        PHARE_WITH_THRUST({ //
            ++sorted;
            thrust::sort(thrust::seq, from, to /*, std::greater<>()*/);
        })
        //
        PHARE_WITH_THRUST_ELSE({ //
            ++sorted;
            std::sort(from, to /*, std::greater<>()*/);
        }) //

        assert(sorted == 1);
    }


    TileSetView<ParticlesTile<Particles, NdArrayView<dim, double>>> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> gap_idx_, add_into_, cap_, left_;

    std::size_t size_;

    Box<int, dim> box_, ghost_box_;
    lobox_t local_ghost_box_;

}; // TileSetSpan


template<typename Particles, std::uint8_t impl_ = 0>
class TileSetVector
{
    using This                    = TileSetVector<Particles, impl_>;
    using SIZE_T                  = unsigned long long int; // cuda issues
    bool static constexpr c_order = true;

    template<typename P, std::uint8_t i>
    friend class TileSetSpan;

    std::uint8_t constexpr static alloc_impl()
    {
        if (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
            return 1;
        return 0;
    }

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr impl_v       = impl_;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr layout_mode  = Particles::layout_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::view_t;
    using per_tile_particles           = Particles;

    template<typename T>
    using vec_helper = PHARE::Vector<T, alloc_mode, alloc_impl()>;

    template<typename T>
    using nd_array_t = NdArrayVector<dim, T, c_order, alloc_mode>;

    using size_t_vector = typename vec_helper<std::size_t>::vector_t;

    using VecTile = ParticlesTile<Particles, nd_array_t<double>>;
    using SpnTile = ParticlesTile<PSpan_t, NdArrayView<dim, double>>;


    TileSetVector(box_t const& box, auto const ghost_cells)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        left_.zero();
        cap_.zero();

        reset_views();
        TileSet<VecTile, alloc_mode>::build_links(particles_);
        TileSet<SpnTile, alloc_mode>::build_links(particles_views_);
    }

    TileSetVector(TileSetVector&& that)
        : ghost_cells_{that.ghost_cells_}
        , box_{that.box_}
        , ghost_box_{that.ghost_box_}
        , particles_{std::move(that.particles_)}
        , total_size{that.total_size}
    {
    }

    TileSetVector(TileSetVector const& that)
        : ghost_cells_{that.ghost_cells_}
        , box_{that.box_}
        , ghost_box_{that.ghost_box_}
        , particles_{that.particles_.copy(TileSetter<dim>{that.box_, that.ghost_cells_})}
        , total_size{that.total_size}
    {
        sync();
    }


    TileSetVector& operator=(TileSetVector&&)      = default;
    TileSetVector& operator=(TileSetVector const&) = default;

    template<auto type = ParticleType::Domain>
    auto& reserve_ppc(std::size_t const& ppc);

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
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(return mkn::gpu::as_super(v));
        else
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



    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    auto& views() { return particles_views_; }
    auto& views() const { return particles_views_; }
    auto& view(locell_t const& cell) { return (*particles_views_.at(cell))(); }

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
    auto local_box(Box<int, dim> const& from) const
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    // PHASE 0 == init
    // PHASE 1 == post move sync
    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync();

    template<auto type>
    void sync_moved();

    template<auto type>
    void trim();



    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        for (auto& tile : particles_views_)
            tile().clear();

        reset_views();
        total_size = 0;
    }




    void static resize(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { resize(v, s, copy); });
        else
            resize(ps.particles_, s, copy);
    }

    template<typename V>
    void static resize(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::resize(v, s, copy));
        else
            v.resize(s);
    }

    void static reserve(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { reserve(v, s, copy); });
        else
            reserve(ps.particles_, s, copy);
    }

    template<typename V>
    void static reserve(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::reserve(v, s, copy));
        else
            v.reserve(s);
    }




    auto local_tile_cell(std::array<int, dim> const& cell) const
    {
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }



protected:
    void sync_check_realloc();

    auto on_tiles(auto&& fn)
    {
        for (auto& tile : particles_)
            fn(tile);
    }


    std::size_t ghost_cells_;

    Box<int, dim> box_, ghost_box_;

    nd_array_t<size_t_vector> gaps_{local_box().shape()};
    nd_array_t<Span<std::size_t>> gap_views_{local_box().shape()};

    nd_array_t<SIZE_T> gap_idx_{local_box().shape()};
    nd_array_t<SIZE_T> add_into_{local_box().shape()};
    nd_array_t<SIZE_T> left_{local_box().shape()};
    nd_array_t<SIZE_T> cap_{local_box().shape()};
    nd_array_t<std::size_t> cell_size_{local_box().shape()};


    // tiles build from amrbox, but `.at()` function maps ghosts box
    TileSet<VecTile, alloc_mode> particles_{TileSetter<dim>{box_, ghost_cells_}, ghost_cells_};
    TileSet<SpnTile, alloc_mode> particles_views_ = generate_from<alloc_mode>(
        TileSetter<dim>{box_, ghost_cells_}, reset_particle_views_fn(), particles_, ghost_cells_);


    std::size_t total_size = 0;

}; // TileSetVector<Particles>




template<typename Particles, std::uint8_t impl>
template<auto type>
auto& TileSetVector<Particles, impl>::reserve_ppc(std::size_t const& ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    std::size_t const additional = ppc < 50 ? 10 : ppc * .2; // 20% overallocate
    std::size_t const buffered   = ppc + additional;

    if constexpr (type == ParticleType::Domain)
        on_tiles([&](auto& tile) {
            reserve(tile(), buffered * tile.size());
            reserve(gaps_(local_cell(tile.lower)), additional * tile.size());
        });

    if constexpr (type == ParticleType::Ghost)
        on_tiles([&](auto& tile) {
            reserve(tile(), additional * tile.size());
            reserve(gaps_(local_cell(tile.lower)), additional * tile.size());
        });

    sync<2>();

    return *this;
}



template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::trim() // change to erase(box)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
}



template<typename Particles, std::uint8_t impl>
void TileSetVector<Particles, impl>::sync_check_realloc()
{
    PHARE_LOG_SCOPE(3, "TileSetVector::sync_check_realloc");

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        auto const bix       = local_cell(particles_[i].lower);
        auto const& nu       = add_into_(bix);
        auto& real           = particles_[i]();
        auto const& left     = gap_idx_(bix);
        auto const& old_size = real.size();
        cell_size_(bix)      = old_size;
        auto const& new_size = real.size() + nu - left;
        reserve(real, real.size() + nu); // we must add new before removing leavers
        resize(real, new_size);
        cap_(bix) = real.capacity();
        assert(real.size() == cell_size_(bix) + nu - left);
        assert(cap_(bix) >= real.size());
        add_into_(bix) = 0;
    };
}


template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::sync_moved()
{
    sync_check_realloc();
    reset_views();

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        auto const bix = local_cell(particles_[i].lower);
        particles_views_[i]().resize(cell_size_(bix));
    };
}

template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void TileSetVector<Particles, impl>::sync()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    PHARE_LOG_SCOPE(3, "TileSetVector::sync");
    // PHARE_LOG_LINE_STR("sync " << magic_enum::enum_name(type));

    if constexpr (PHASE == 1)
        sync_check_realloc();

    total_size = 0;
    for (auto const& tile : particles_)
    {
        total_size += tile().size();
        cell_size_(local_cell(tile.lower)) = tile().size();
    }

    auto const per_cell = [&](auto& tile) {
        auto const bix  = local_cell(tile.lower);
        auto const& cs  = cell_size_(bix);
        auto const& cap = tile().capacity();
        auto& gaps      = gaps_(bix);
        if (gaps.size() < cs)
        {
            reserve(gaps, cap, false);
            resize(gaps, cs, false);
        }
        cap_(bix) = cap;
    };

    on_tiles(per_cell);

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

    template<std::size_t size>
    using array_type = per_tile_particles::template array_type<size>;

    auto static constexpr impl_v       = Super::impl_v;
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
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    template<typename T>
    struct iterator_impl;

    template<typename T, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto erase(iterator_impl<T> a, iterator_impl<T> b)
    {
        return Super::erase(a, b);
    }

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

    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }


    auto& icell_changer(std::array<std::uint32_t, dimension> const& cell, std::size_t const idx,
                        std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        bool constexpr static ATOMIC = true;
        bool constexpr static GPU    = alloc_mode == AllocatorMode::GPU_UNIFIED;
        using Op                     = Operators<typename Super::SIZE_T, ATOMIC, GPU>;

        auto& gidx      = Super::gap_idx_(cell);
        auto const nidx = Op{gidx}.increment_return_old();
        auto& gaps      = Super::gaps_(cell);
        assert(nidx < gaps.size());
        gaps[nidx] = idx;

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
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return index_wrapper<This const>{this, s};
    }

    void print() const {}
    void check() const
    {
        // for (auto const& tile : Super::particles_())
        // {
        //     auto const gb = tile.field_box();

        //     for (auto const& particle : tile())
        //     {
        //         assert(isIn(particle, gb));
        //     }
        // }

        // if constexpr (storage_mode == StorageMode::VECTOR)
        //     for (auto const& tile : Super::particles_views_())
        //     {
        //         auto const gb = tile.field_box();

        //         for (auto const& particle : tile())
        //         {
        //             assert(isIn(particle, gb));
        //         }
        //     }
    }

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
        throw std::runtime_error("finish this");

        return 10;
    }


}; // TileSetParticles<Super>


template<typename Particles>
struct tile_set_iterator_base
{
    tile_set_iterator_base(Particles ps)
        : particles{ps}
    {
    }

    Particles particles;

    std::size_t l = 0, i = 0;
};

template<auto layout_mode, typename Particles>
struct tile_set_iterator_storage;

#if PHARE_HAVE_THRUST
template<typename Particles>
struct tile_set_iterator_storage<LayoutMode::SoA, Particles>
    : public tile_set_iterator_base<Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using Super                    = tile_set_iterator_base<Particles>;
    using per_tile_particles       = typename std::decay_t<Particles>::per_tile_particles;
    using Particle_t = typename SoAZipParticle_t<per_tile_particles, is_const>::value_type;
    using Super::l, Super::i, Super::particles;

    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    void set() { particle = Particle_t{particles.data()[l], i}; }

    Particle_t particle;
};

template<typename Particles>
struct tile_set_iterator_storage<LayoutMode::SoAVX, Particles>
    : public tile_set_iterator_base<Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using Super                    = tile_set_iterator_base<Particles>;
    using per_tile_particles       = std::decay_t<Particles>::per_tile_particles;
    using Particle_t = typename SoAVXZipParticle_t<per_tile_particles, is_const>::value_type;
    using Super::l, Super::i, Super::particles;

    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    void set() { particle = Particle_t{particles.data()[l], i}; }

    Particle_t particle;
};

#else

#endif // PHARE_HAVE_THRUST

template<typename Particles>
struct tile_set_iterator_storage<LayoutMode::AoS, Particles>
    : public tile_set_iterator_base<Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using Super                    = tile_set_iterator_base<Particles>;
    using per_tile_particles       = std::decay_t<Particles>::per_tile_particles;
    using Particle_t               = per_tile_particles::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;
    using Super::l, Super::i, Super::particles;

    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    auto& charge() _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    auto& charge() const _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() const _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() const _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() const _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() const _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    void set() { /*noop*/ }
};

template<auto layout_mode, typename Particles>
struct tile_set_iterator_storage : public tile_set_iterator_base<Particles>
{
    using Super = tile_set_iterator_base<Particles>;

    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps}
    {
        throw std::runtime_error("no");
    }
};

template<typename T>
struct tile_set_iterator_super
{
    using per_tile_particles = typename std::decay_t<T>::per_tile_particles;
    using value_type         = tile_set_iterator_storage<per_tile_particles::layout_mode, T>;
};
template<typename T>
using tile_set_iterator_super_v = typename tile_set_iterator_super<T>::value_type;

template<typename OuterSuper>
template<typename T>
struct TileSetParticles<OuterSuper>::iterator_impl : public tile_set_iterator_super_v<T>
{
    auto static constexpr dimension = OuterSuper::dimension;
    using Super                     = tile_set_iterator_super_v<T>;
    using outer_type                = std::decay_t<T>;
    using difference_type           = std::size_t;
    using iterator_category         = std::forward_iterator_tag;
    using Particle_t                = typename OuterSuper::Particle_t;
    using value_type                = Particle_t;
    using pointer                   = Particle_t*;
    using reference                 = Particle_t&;
    using Super::l, Super::i, Super::particles;

    iterator_impl(T& particles_, bool end = false)
        : Super{particles_}
    {
        if (end)
        {
            l = particles.particles_.size() - 1;
            i = particles.size(l);
        }
        else
            inc(); // find first particle
    }
    iterator_impl(iterator_impl&& that)                 = default;
    iterator_impl(iterator_impl const& that)            = default;
    iterator_impl& operator=(iterator_impl&& that)      = default;
    iterator_impl& operator=(iterator_impl const& that) = default;

    auto end() { return iterator_impl{particles, true}; }

    auto& set()
    {
        Super::set();
        return *this;
    }
    auto& inc()
    {
        i = 0;
        for (; l < particles.particles_.size(); ++l)
            if (*this == this->end() or particles.size(l) > 0)
                break;
        return set();
    }

    auto& operator++()
    {
        auto const last = end();
        if (*this == last)
            return *this;
        ++i;
        if (*this == last)
            return *this;
        if (i == particles.size(l))
        {
            ++l;
            return inc();
        }
        return set();
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto operator==(iterator_impl const& that) const { return l == that.l and i == that.i; }
    auto operator!=(iterator_impl const& that) const { return !(*this == that); }
    auto operator<(iterator_impl const& that) const
    {
        // PHARE_LOG_LINE_STR(l << " " << i);
        if (l < that.l)
            return true;
        if (l == that.l)
            return i < that.i;
        return false;
    }


    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }
};




template<typename Particles, std::uint8_t impl>
template<typename... Args>
void TileSetSpan<Particles, impl>::sync(Args&&... args) _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(3, "TileSetSpan::sync(stream)");

    auto& view = *this;

    if constexpr (Particles::alloc_mode == AllocatorMode::CPU)
    {
        for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
            sync_tile_add_new(tidx);

        for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
            sync_tile_rm_left(tidx);
    }
    else if (Particles::alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        if (!PHARE_HAVE_MKN_GPU)
            throw std::runtime_error("no gpu impl");

        PHARE_WITH_MKN_GPU({
            auto const& [stream] = std::forward_as_tuple(args...);
            mkn::gpu::GDLauncher<true>{particles_.size()} //
                .stream(stream, [=] _PHARE_ALL_FN_() mutable {
                    view.sync_tile_add_new(mkn::gpu::idx());
                    __syncthreads();
                    view.sync_tile_rm_left(mkn::gpu::idx());
                });
        })
    }
    else
        throw std::runtime_error("no");
}


template<typename Particles, std::uint8_t impl>
void TileSetSpan<Particles, impl>::sync_tile_add_new(std::size_t const tidx) _PHARE_ALL_FN_
{
    bool constexpr static ATOMIC = true;
    bool constexpr static GPU    = alloc_mode == AllocatorMode::GPU_UNIFIED;

    using Op = Operators<SIZE_T, ATOMIC, GPU>;

    auto& tile     = particles_[tidx];
    auto& real     = tile();
    auto const bix = local_cell(tile.lower);

    auto const& n_gaps = gap_idx_(bix);
    {
        auto& gaps = gaps_(bix);
        sort(gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
    }

    auto& left       = left_(bix);
    auto const& gaps = gaps_(bix);
    for (std::size_t i = 0; i < n_gaps; ++i)
    {
        auto const& gidx    = gaps[n_gaps - (1 + i)];
        auto const& newcell = local_tile_cell(real.iCell(gidx));
        auto const& cap     = cap_(newcell);
        auto& ntile         = (*particles_.at(newcell));
        auto& nparts        = ntile();
        auto npidx          = nparts.size();

        while (true)
        {
            assert(npidx < cap);

            auto inc = npidx + 1;
            auto old = Op::compare_and_swap(nparts.size_address(), npidx, inc);
            if (npidx != old)
            {
                ++npidx;
                continue;
            }
            else
                break;
        }

        assert(not isIn(real[gidx], tile));
        assert(isIn(real[gidx], ntile));
        nparts.assign(real, gidx, npidx);
        assert(isIn(nparts[npidx], ntile));
        ++left;
    }

    assert(left == n_gaps);
}


template<typename Particles, std::uint8_t impl>
void TileSetSpan<Particles, impl>::sync_tile_rm_left(std::size_t const tidx) _PHARE_ALL_FN_
{
    auto& tile       = particles_[tidx];
    auto& real       = tile();
    auto const bix   = local_cell(tile.lower);
    auto const& gaps = gaps_(bix);
    auto& left       = left_(bix);
    auto& gaps_size  = gap_idx_(bix);

    while (left)
    {
        auto const& rsize = real.size();
        auto const& pidx  = gaps[gaps_size - 1];

        PHARE_DEBUG_DO({
            auto const& p = real[pidx];
            assert(not isIn(p, tile));
        })

        if (pidx != rsize - 1)
            real.assign(rsize - 1, pidx);

        real.pop_back();
        --gaps_size;
        --left;
    }

    PHARE_ASSERT(left == 0);
    PHARE_ASSERT(gaps_size == 0);
}


template<auto type>
void sync_aos_ts(auto& particles, auto&&... args)
{
    particles.template sync_moved<type>(); // realloc
    (*particles).sync(args...);            // swap

    // assert final size
    PHARE_DEBUG_DO({
        auto& tiles = particles();
        auto& views = particles.views();

        for (std::size_t i = 0; i < particles().size(); ++i)
        {
            auto const& real_size = tiles[i]().size();
            auto const& view_size = views[i]().size();

            // see
            assert(view_size == real_size);
        }
    })

    particles.template sync<0, type>(); // finalize
}




} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP */
