#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP


#include "core/operators.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"

#include "core/utilities/span.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "initializer/data_provider.hpp"


#include <tuple>
#include <iterator>
#include <stdexcept>


namespace PHARE::core
{

struct TileSetParticleArrayDetails : public ParticleArrayDetails
{
    std::size_t const interp_order = 1; // for fields per tile
    std::size_t const tile_size    = 4;

    template<typename GridLayout_t>
    TileSetParticleArrayDetails static FROM(GridLayout_t const& layout,
                                            initializer::PHAREDict const& dict)
    {
        auto const super = ParticleArrayDetails::FROM(layout);
        TileSetParticleArrayDetails const defaults{};
        return {
            {super},
            cppdict::get_value(dict, "interp_order", defaults.interp_order),
            cppdict::get_value(dict, "tile_size", defaults.tile_size),
        };
    }
};




template<typename Particles, typename NdArray_t>
class ParticlesTile : public Box<std::int32_t, Particles::dimension>
{
    auto constexpr static dim = Particles::dimension;
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

    template<typename... Args>
    auto constexpr static bounded_ghost_box(Args&&... args)
    {
        auto const& [tilebox, particle_array] = std::forward_as_tuple(args...);
        auto tgbox                            = grow(tilebox, 2); // assume interp 1 for now
        auto bounded_tgbox                    = tgbox * particle_array.safe_box();
        assert(bounded_tgbox);
        (*bounded_tgbox).upper += 1; // +1 for primal
        return *bounded_tgbox;
    }

public:
    template<typename Ps, typename Nd>
    ParticlesTile(ParticlesTile<Ps, Nd>& tile) // span constructor
        : Super{tile}
        , particles{tile(), tile().size()}
        , field_ghost_box{tile.field_ghost_box}
        , rho{*tile.rho}
        , fx{*tile.fx}
        , fy{*tile.fy}
        , fz{*tile.fz}
    {
    }

    template<typename... Args>
    ParticlesTile(Super const& box, Args&&... args)
        : Super{box}
        , particles{[&]() {
            if constexpr (std::is_constructible_v<Particles, Super>)
                return Particles{box};
            else
                return Particles{};
        }()}
        , field_ghost_box{bounded_ghost_box(box, args...)}
        , rho{ndarray_builder(field_ghost_box)}
        , fx{ndarray_builder(field_ghost_box)}
        , fy{ndarray_builder(field_ghost_box)}
        , fz{ndarray_builder(field_ghost_box)}
    {
    }


    ParticlesTile(ParticlesTile<Particles, NdArray_t>&) = default;


    ParticlesTile(ParticlesTile&&)                 = default;
    ParticlesTile(ParticlesTile const&)            = default;
    ParticlesTile& operator=(ParticlesTile&&)      = default;
    ParticlesTile& operator=(ParticlesTile const&) = default;

    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

    auto& operator()() _PHARE_ALL_FN_ { return particles; }
    auto& operator()() const _PHARE_ALL_FN_ { return particles; }

    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }

    auto fields() _PHARE_ALL_FN_ { return std::forward_as_tuple(rho, fx, fy, fz); }
    auto fields() const _PHARE_ALL_FN_ { return std::forward_as_tuple(rho, fx, fy, fz); }
    auto& field_box() const _PHARE_ALL_FN_ { return field_ghost_box; }

    template<typename P, typename N>
    void reset(ParticlesTile<P, N>& tile)
    {
        (*this)().reset(tile());
        // fx.reset(tile.fx);
    }

private:
    Particles particles;
    Super field_ghost_box; // tilebox + field ghosts
    NdArray_t rho, fx, fy, fz;

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
        {
            arr.check();
            return arr.particles_views_.make_view();
        }
    }

    template<typename TileSetArray>
    auto resolve_gaps(TileSetArray& arr)
    {
        if constexpr (TileSetArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
        {
            arr.check();
            return *arr.gap_views_;
        }
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;

    template<typename TileSetArray>
    TileSetSpan(TileSetArray& arr)
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , off_sets_{arr.off_sets_}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cap_{arr.cap_}
        , left_{arr.left_}
        , p2c_{arr.p2c_.data(), arr.p2c_.size()}
        , size_{arr.size()}
        , box_{arr.box_}
        , ghost_box_{arr.ghost_box()}
        , safe_box_{arr.safe_box_}
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
    auto& safe_box() const { return safe_box_; }
    auto local_cell(std::array<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return as_local_cell(safe_box_, icell);
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

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain, typename... Args>
    void sync(Args&&... args) _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_add_new() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left() _PHARE_ALL_FN_;

    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        size_ = 0;
    }

    void reset_p2c(std::array<std::uint32_t, dim>* cells, std::size_t size)
    {
        p2c_ = Span<std::array<std::uint32_t, dim>>{cells, size};
    }

    auto local_tile_cell(std::array<int, dim> const& cell) const _PHARE_ALL_FN_
    {
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

protected:
    TileSetView<ParticlesTile<Particles, NdArrayView<dim, double>>> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> off_sets_, gap_idx_, add_into_, cap_, left_;
    Span<std::array<std::uint32_t, dim>> p2c_;
    std::size_t size_;

    Box<int, dim> box_, ghost_box_, safe_box_;
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
    using vec_helper = PHARE::Vector<T, alloc_mode, 1>;

    template<typename T>
    using nd_array_t = NdArrayVector<dim, T, c_order, alloc_mode>;

    using size_t_vector = typename vec_helper<std::size_t>::vector_t;

    using VecTile = ParticlesTile<Particles, nd_array_t<double>>;
    using SpnTile = ParticlesTile<PSpan_t, NdArrayView<dim, double>>;


    TileSetVector(box_t const& box /*= {}*/, TileSetParticleArrayDetails const& deets = {})
        : details{deets}
        , ghost_cells_{deets.ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, deets.ghost_cells)}
        , safe_box_{grow(box, deets.ghost_cells + 1)}
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

    TileSetVector(TileSetVector const& from)            = default;
    TileSetVector(TileSetVector&& from)                 = default;
    TileSetVector& operator=(TileSetVector&& from)      = default;
    TileSetVector& operator=(TileSetVector const& from) = default;

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
    void /*Particle_t&*/ emplace_back(Particle_t const& p)
    {
        // auto& np = get_vec(particles_(local_cell(p.iCell()))).emplace_back(p);
        auto const locell = local_cell(p.iCell());
        /*auto& np = */ (*particles_.at(local_cell(p.iCell())))().emplace_back(p);
        if constexpr (inc_)
            _inc(locell);
        // return np;
    }

    template<bool inc_ = true>
    void emplace_back(Particles& dst, Particles const& src, std::size_t const& idx)
    {
        dst.emplace_back(src, idx);
    }

    template<typename V>
    static auto& get_vec(V& v)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(return mkn::gpu::as_super(v));
        }
        else
        {
            return v;
        }
    }

    void push_back(Particle_t&& p) { emplace_back(p); }
    void push_back(Particle_t const& p) { emplace_back(p); }


    void reset_views()
    {
        // update_from(reset_particle_views_fn(), particles_views_);
        update_from(reset_gap_views_fn(), gap_views_);
        for (std::size_t i = 0; i < particles_.size(); ++i)
            particles_views_[i]().reset(particles_[i]());
    }
    auto reset_gap_views_fn()
    {
        return [&](auto const i) { return make_span(*(gaps_.data() + i)); };
    }
    auto reset_particle_views_fn()
    {
        return [&](std::size_t const i) { return SpnTile{particles_[i]}; };
    }

    void reset_index_wrapper_map();

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }
    auto& safe_box() const { return safe_box_; } // allocating box

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return (*particles_.at(cell))(); }
    // auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    auto& views() { return particles_views_; }
    auto& view(locell_t const& cell) { return (*particles_views_.at(cell))(); }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(safe_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }
    auto local_box() const
    {
        return box_from_zero_to_upper_minus_one(
            safe_box_.shape().template toArray<std::uint32_t>());
    }
    auto local_box(Box<int, dim> const& from) const
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync(); // after move

    template<auto type>
    void trim();



    void clear()
    {
        for (auto& tile : particles_)
            tile().clear();
        for (auto& tile : particles_views_)
            tile().clear();
        total_size = 0;
    }


    auto& insert(TileSetVector const& src);

    template<typename GridLayout_t>
    auto& insert_domain_from(TileSetVector const& src, GridLayout_t const& layout);


    void replace_from(This const& that)
    {
        throw std::runtime_error("fix");
        if (this == &that) // just in case
            return;

        cell_size_ = that.cell_size_;
        particles_ = that.particles_;
        box_       = that.box_;
        ghost_box_ = that.ghost_box_;
        // p2c_       = that.p2c_;
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

    template<typename View_t>
    void reset_p2c(View_t& view)
    {
        view.reset_p2c(p2c_.data(), p2c_.size());
    }


    auto local_tile_cell(std::array<int, dim> const& cell) const
    {
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

protected:
    TileSetParticleArrayDetails const details;

    std::size_t ghost_cells_;

    Box<int, dim> box_, ghost_box_, safe_box_;

    typename vec_helper<std::array<std::uint32_t, dim>>::vector_t p2c_;
    nd_array_t<SIZE_T> off_sets_{local_box().shape()};

    nd_array_t<size_t_vector> gaps_{local_box().shape()};
    nd_array_t<Span<std::size_t>> gap_views_{local_box().shape()};

    nd_array_t<SIZE_T> gap_idx_{local_box().shape()};
    nd_array_t<SIZE_T> add_into_{local_box().shape()};
    nd_array_t<SIZE_T> left_{local_box().shape()};
    nd_array_t<SIZE_T> cap_{local_box().shape()};
    nd_array_t<std::size_t> cell_size_{local_box().shape()};

    TileSet<VecTile, alloc_mode> particles_{safe_box_, details.tile_size, *this};

    // only used for GPU
    TileSet<SpnTile, alloc_mode> particles_views_{
        generate_from<alloc_mode>(reset_particle_views_fn(), particles_, *this)};


    std::size_t total_size = 0;

    template<auto type>
    void sync_cpu_gaps_and_tmp();
    template<auto type>
    void sync_gpu_gaps_and_tmp();

    template<auto type>
    void sync_gpu_gaps_and_tmp_impl0();

    template<auto type>
    void sync_gpu_gaps_and_tmp_impl1();



    auto on_tiles(auto&& fn)
    {
        for (auto& tile : particles_)
            fn(tile);
    }



}; // TileSetVector<Particles>



template<typename Particles, std::uint8_t impl>
auto& TileSetVector<Particles, impl>::insert(TileSetVector const& /*src*/)
{
    // std::size_t added = 0;
    // for (auto const& bix : local_box(box()))
    // {
    //     auto& from = src(bix);
    //     auto& to   = (*this)(bix);
    //     added += from.size();
    //     to.reserve(to.size() + from.size());
    //     std::copy(from.begin(), from.end(), std::back_inserter(to));
    // }
    // if (added)
    //     sync<2>();

    return *this;
}

template<typename Particles, std::uint8_t impl>
template<typename GridLayout_t>
auto& TileSetVector<Particles, impl>::insert_domain_from(TileSetVector const& src,
                                                         GridLayout_t const& layout)
{
    using Partitioner = ParticleArrayPartitioner<PSpan_t>;
    std::size_t added = 1;

    auto const per_tile = [&](auto src_tile) {
        //
        auto const in_domain = Partitioner{src_tile()}(box_);
        auto remaining       = in_domain;

        if (remaining.size() == 0)
            return;

        auto const per_neighbour = [&](auto neigh) {
            if (remaining.size() == 0)
                return;

            auto const not_in_neighbour = Partitioner{src_tile()}.notIn(neigh);
            auto const start            = not_in_neighbour.size();
            auto const end              = remaining.size();

            if (auto const size = end - start)
            {
                auto& dst_tile = *this->particles_.at(neigh.lower - src.safe_box().lower);
                append_particles<ParticleType::Domain>(src_tile().view(start, size), dst_tile(),
                                                       layout);
            }

            remaining = not_in_neighbour;
            src_tile().resize(remaining.size());
        };

        per_neighbour(src_tile);
        traverse_tile_neighbours(src.particles_views_, src_tile, per_neighbour);
    };

    traverse_ghost_boundary_tiles(src.particles_views_, per_tile);

    if (added)
        sync<2>();
    return *this;
}

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
    {
        on_tiles([&](auto& tile) {
            reserve(tile(), additional * tile.size());
            reserve(gaps_(local_cell(tile.lower)), additional * tile.size());
        });
    }

    return *this;
}



template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::trim()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    // currently using views for gpu partition
    reset_views();
    // update_from(reset_particle_views_fn(), particles_views_);

    if constexpr (type == ParticleType::Ghost)
    {
        for (std::size_t i = 0; i < particles_.size(); i++)
        {
            if (particles_views_[i]().size() == 0)
                continue;

            auto& tile = particles_[i];
            auto& tbox = *tile;
            auto& view = particles_views_[i];
            if (tbox * box_ != tbox) // includes ghost cells
                resize(tile(), ParticleArrayPartitioner<PSpan_t>{view()}(box_).size());
        }

        // for (auto const& ghost_layer_box : local_box().remove(local_box(box_)))
        //     for (auto const& bix : ghost_layer_box)
        //     {
        //         cell_size_(bix) = 0;
        //         particles_(bix).clear();
        //         gaps_(bix).clear();
        //     }
    }

    if constexpr (type == ParticleType::Domain)
    {
        throw std::runtime_error("NO");
    }
}

template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::sync_cpu_gaps_and_tmp()
{
    // PHARE_LOG_LINE_STR("sync_cpu_gaps_and_tmp " << magic_enum::enum_name(type));

    for (auto& tile : particles_)
    {
        auto const bix = local_cell(tile.lower);
        auto& real     = tile();
        auto& gaps     = gaps_(bix);
        while (gaps.size())
        {
            real.assign(real.size() - 1, gaps.back());
            gaps.pop_back();
            real.pop_back();
        }
        PHARE_ASSERT(gaps.size() == 0);
    }
}


template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::sync_gpu_gaps_and_tmp_impl0()
{
    { // calculate reserve size
        PHARE_LOG_SCOPE(1, "TileSetVector::sync_gpu_gaps_and_tmp::reserve_scan ");

        for (std::size_t i = 0; i < particles_.size(); ++i)
        {
            auto const bix = local_cell(particles_[i].lower);
            auto& real     = particles_[i]();
            if constexpr (any_in(impl_v, 1, 2))
            {
                // push_back is done on device requires over allocation
                // cell_size_(bix) = particles_views_(bix).size();
                // real.resize(particles_views_(bix).size());
                resize(real, particles_views_[i]().size());
            }
            reserve(real, real.size() + add_into_(bix));
        };
    }

    { // add incoming particles
        PHARE_LOG_SCOPE(1, "TileSetVector::sync_gpu_gaps_and_tmp::add ");

        for (std::size_t i = 0; i < particles_.size(); ++i)
        {
            auto const bix        = local_cell(particles_[i].lower);
            auto& real            = particles_[i]();
            auto const& gaps      = gaps_(bix);
            auto const& gaps_size = gap_idx_(bix);
            for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
            {
                auto const& idx     = gaps[gidx];
                auto const& newcell = local_cell(real.iCell(idx));
                emplace_back((*particles_.at(newcell))(), real, idx);
            }
        }
    }

    { // delete outgoing particles
        PHARE_LOG_SCOPE(1, "TileSetVector::sync_gpu_gaps_and_tmp::delete ");
        for (std::size_t i = 0; i < particles_.size(); ++i)
        {
            auto const bix        = local_cell(particles_[i].lower);
            auto const& gaps_size = gap_idx_(bix);
            {
                auto& gaps = gaps_(bix);
                std::sort(gaps.begin(), gaps.begin() + gaps_size, std::greater<>()); // use thrust ?
            }
            auto& real       = particles_[i]();
            auto const& gaps = gaps_(bix);
            for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
            {
                auto const& idx = gaps[gidx];
                real.assign(real.size() - 1, idx);
                real.pop_back();
            }
            gap_idx_(bix)  = 0;
            add_into_(bix) = 0;
        }
    }
}

template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::sync_gpu_gaps_and_tmp_impl1()
{
}



template<typename Particles, std::uint8_t impl>
template<auto type>
void TileSetVector<Particles, impl>::sync_gpu_gaps_and_tmp()
{
    // PHARE_LOG_LINE_STR("sync_gpu_gaps_and_tmp " << magic_enum::enum_name(type));
    PHARE_LOG_SCOPE(1, "TileSetVector::sync_gpu_gaps_and_tmp ");

    // if constexpr (impl == 0)
    // {
    sync_gpu_gaps_and_tmp_impl0<type>();
    // }
    // else if constexpr (impl == 1)
    // {
    //     sync_gpu_gaps_and_tmp_impl1<type>();
    // }
    // else
    //     throw std::runtime_error("No impl");
}

template<typename Particles, std::uint8_t impl>
void TileSetVector<Particles, impl>::reset_index_wrapper_map()
{
    resize(p2c_, total_size);

    auto const fill = [](auto p, auto o, auto s, auto b) { std::fill(p + o, p + o + s, b); };

    std::size_t offset = 0;

    on_tiles([&](auto& tile) {
        auto const bix = local_cell(tile.lower);

        auto const& cs = cell_size_(bix);
        resize(gaps_(bix), cs);
        off_sets_(bix) = offset;

        if (cs)
        {
            if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
            {
                PHARE_WITH_THRUST( //
                    thrust::fill(thrust::device, p2c_.begin() + offset, p2c_.begin() + offset + cs,
                                 *bix));
                PHARE_WITH_THRUST_ELSE(
                    PHARE_LOG_LINE_SS("Thrust not found for TileSetVector<Particles, "
                                      "impl>::reset_index_wrapper_map"); //
                    fill(p2c_.begin(), offset, cs, bix);                 //
                )
            }
            else
                fill(p2c_.begin(), offset, cs, bix);
        }

        offset += cs;
        cap_(bix) = tile().capacity();
    });
}

template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void TileSetVector<Particles, impl>::sync()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    PHARE_LOG_SCOPE(1, "TileSetVector::sync");
    // PHARE_LOG_LINE_STR("sync " << static_cast<std::uint32_t>(PHASE) << " "
    //                            << magic_enum::enum_name(type));
    // auto const lbox = local_box();

    if constexpr (PHASE < 2 and alloc_mode == AllocatorMode::CPU)
        sync_cpu_gaps_and_tmp<type>();

    if constexpr (PHASE < 2 and alloc_mode == AllocatorMode::GPU_UNIFIED)
        sync_gpu_gaps_and_tmp<type>();

    if constexpr (PHASE == 1 || type == ParticleType::Domain)
        trim<ParticleType::Ghost>();

    total_size = 0;
    for (auto const& tile : particles_)
    {
        total_size += tile().size();
        cell_size_(local_cell(tile.lower)) = tile().size();
    }

    if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        PHARE_LOG_SCOPE(1, "TileSetVector::sync::reset");
        static_assert(impl < 3); // otherwise unhandled
        if constexpr (impl < 2)
        {
            // reset_index_wrapper_map();
        }
        else // if (PHASE == 2)
        {
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

            if constexpr (type == ParticleType::Domain)
                on_tiles(per_cell);

            // if constexpr (ParticleType::Domain)
            // {}
        }
    }

    {
        PHARE_LOG_SCOPE(1, "TileSetVector::sync::reset_views");
        reset_views();
    }
}


template<typename Super_>
struct TileSetParticles : public Super_
{
    using Super              = Super_;
    using This               = TileSetParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_tile_particles = typename Super::per_tile_particles;

    auto static constexpr impl_v       = Super::impl_v;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::local_box;
    using Super::particles_;
    using Super::size;


    template<typename... Args>
    TileSetParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }



    TileSetParticles(TileSetParticles const& from)            = default;
    TileSetParticles(TileSetParticles&& from)                 = default;
    TileSetParticles& operator=(TileSetParticles&& from)      = default;
    TileSetParticles& operator=(TileSetParticles const& from) = default;

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



    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<std::uint32_t, dimension> const& local) const
    {
        std::vector<Particle_t> out;
        out.reserve(sum_from(local, [&](auto const& b) { return particles_(b.toArray()).size(); }));
        for (auto const& b : local)
            std::copy(particles_(b).begin(), particles_(b).end(), std::back_inserter(out));
        return out;
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<int, dimension> const& amr) const
    {
        return (*this)[Super::local_box(amr)];
    }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }


    template<typename _Particle>
    auto& icell_changer(_Particle const& p, std::array<std::uint32_t, dimension> const& cell,
                        std::size_t const idx,
                        std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            Super::get_vec(Super::gaps_(cell)).emplace_back(idx);
            if (isIn(newcell, Super::ghost_box()))
                Super::get_vec((*this)(Super::local_tile_cell(newcell))).emplace_back(p).iCell()
                    = newcell;
        }
        else if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            using Op = Operators<typename Super::SIZE_T, true>;

            Super::gaps_(cell)[Op{Super::gap_idx_(cell)}.increment_return_old()] = idx;

            if (isIn(newcell, Super::ghost_box()))
            {
                auto const nc = Super::local_tile_cell(newcell);
                Op{Super::add_into_(nc)}.increment_return_old();
            }
        }
        else
            throw std::runtime_error("no");

        return *this;
    }


    template<typename T>
    struct index_wrapper;

    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return index_wrapper<This const>{this, s};
    }

    void print() const {}

    void check() const {}

    auto max_size() const
    {
        return max_from(this->particles_,
                        [](auto const& v, auto const& i) { return v.data()[i].size(); });
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

    // auto& operator*() _PHARE_ALL_FN_ { return particle; }
    // auto& operator*() const _PHARE_ALL_FN_ { return particle; }

    // per_tile_particles* particles;
    void set() { particle = Particle_t{particles.data()[l], i}; }

    Particle_t particle;
};

template<typename Particles>
struct tile_set_iterator_storage<LayoutMode::SoAVX, Particles>
    : public tile_set_iterator_base<Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using Super                    = tile_set_iterator_base<Particles>;
    using per_tile_particles       = typename std::decay_t<Particles>::per_tile_particles;
    using Particle_t = typename SoAVXZipParticle_t<per_tile_particles, is_const>::value_type;
    using Super::l, Super::i, Super::particles;


    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    // auto& operator*() _PHARE_ALL_FN_ { return particle; }
    // auto& operator*() const _PHARE_ALL_FN_ { return particle; }

    void set() { particle = Particle_t{particles.data()[l], i}; }

    // auto& charge() const _PHARE_ALL_FN_ { return particle.charge(); }
    // auto& delta() const _PHARE_ALL_FN_ { return particle.delta(); }
    // auto& iCell() const _PHARE_ALL_FN_ { return particle.iCell(); }
    // auto& weight() const _PHARE_ALL_FN_ { return particle.weight(); }
    // auto& v() const _PHARE_ALL_FN_ { return particle.v(); }

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
    using per_tile_particles       = typename std::decay_t<Particles>::per_tile_particles;
    using Particle_t               = typename Particles::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;
    using Super::l, Super::i, Super::particles;

    tile_set_iterator_storage(Particles ps) _PHARE_ALL_FN_ : Super{ps} {}

    // auto& operator*() _PHARE_ALL_FN_ { return *particle; }
    // auto& operator*() const _PHARE_ALL_FN_ { return *particle; }

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

    // Particle_p particle;
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



    // auto& operator*() _PHARE_ALL_FN_ { return particles.data()[l][i]; }
    // auto& operator*() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }


    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    // T particles;
    // std::size_t l = 0, i = 0;
};


template<auto layout_mode, typename Particles>
struct tile_set_index_wrapper_storage;

#if PHARE_HAVE_THRUST
template<typename Particles>
struct tile_set_index_wrapper_storage<LayoutMode::SoA, Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using per_tile_particles       = typename std::decay_t<Particles>::per_tile_particles;
    using Particle_t = typename SoAZipParticle_t<per_tile_particles, is_const>::value_type;

    template<typename PerTileParticles_t>
    tile_set_index_wrapper_storage(PerTileParticles_t p, std::size_t const i) _PHARE_ALL_FN_
        : /*particles{p},*/
          particle{*p, i}
    {
    }

    auto& operator*() _PHARE_ALL_FN_ { return particle; }
    auto& operator*() const _PHARE_ALL_FN_ { return particle; }

    per_tile_particles* particles;
    Particle_t particle;
};
#else

#endif // PHARE_HAVE_THRUST

template<typename Particles>
struct tile_set_index_wrapper_storage<LayoutMode::AoS, Particles>
{
    bool static constexpr is_const = std::is_const_v<std::remove_reference_t<Particles>>;
    using per_tile_particles       = typename std::decay_t<Particles>::per_tile_particles;
    using Particle_t               = typename Particles::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;

    template<typename PerTileParticles_t>
    tile_set_index_wrapper_storage(PerTileParticles_t p, std::size_t const i) _PHARE_ALL_FN_
        : /*particles{p},*/
          particle{&p->data()[i]}
    {
    }

    auto& operator*() _PHARE_ALL_FN_ { return *particle; }
    auto& operator*() const _PHARE_ALL_FN_ { return *particle; }

    // per_tile_particles* particles;
    Particle_p particle;
};

template<auto layout_mode, typename Particles>
struct tile_set_index_wrapper_storage
{
    // unused
}; // default;

template<typename T>
struct tile_set_index_wrapper_super
{
    using per_tile_particles = typename std::decay_t<T>::per_tile_particles;
    using value_type         = tile_set_index_wrapper_storage<per_tile_particles::layout_mode, T>;
};


template<typename ParticlesSuper>
template<typename T>
struct TileSetParticles<ParticlesSuper>::index_wrapper
    : public tile_set_index_wrapper_super<T>::value_type
{
    using outer_t = std::decay_t<T>;
    using Super   = typename tile_set_index_wrapper_super<T>::value_type;

    auto static constexpr dimension = ParticlesSuper::dimension;
    bool static constexpr is_const  = std::is_const_v<std::remove_reference_t<T>>;
    using Particle_t                = typename outer_t::Particle_t;
    using Particle_p = std::conditional_t<is_const, Particle_t const* const, Particle_t*>;


    index_wrapper(T* ts_particles, std::size_t idx_) _PHARE_ALL_FN_
        : Super{&(*ts_particles)(cell(ts_particles, idx_)), index(ts_particles, idx_)},
          ts_particles_ptr{ts_particles},
          idx{idx_}
    {
        PHARE_ASSERT((**this).iCell()[0] > -10 and (**this).iCell()[0] < 1000); // bad memory
        if constexpr (dimension > 1)
            PHARE_ASSERT((**this).iCell()[1] > -10 and (**this).iCell()[1] < 1000); // bad memory
        if constexpr (dimension > 2)
            PHARE_ASSERT((**this).iCell()[2] > -10 and (**this).iCell()[2] < 1000); // bad memory
    }

    auto& c() const _PHARE_ALL_FN_ { return cell(ts_particles_ptr, idx); }
    static auto& cell(T* arr, std::size_t const idx) _PHARE_ALL_FN_ { return arr->p2c_[idx]; }
    auto i() const _PHARE_ALL_FN_ { return index(ts_particles_ptr, idx); }
    static auto index(T* arr, std::size_t const idx) _PHARE_ALL_FN_
    {
        return idx - arr->off_sets_(cell(arr, idx));
    }


    auto icell_changer(std::array<int, dimension> const& newcell) _PHARE_ALL_FN_
    {
        ts_particles_ptr->icell_changer(**this, c(), i(), newcell);
    }


    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }

    auto& operator*() _PHARE_ALL_FN_ { return *super(); }
    auto& operator*() const _PHARE_ALL_FN_ { return *super(); }

    auto& charge() const _PHARE_ALL_FN_ { return (**this).charge(); }
    auto& delta() const _PHARE_ALL_FN_ { return (**this).delta(); }
    auto& iCell() const _PHARE_ALL_FN_ { return (**this).iCell(); }
    auto& weight() const _PHARE_ALL_FN_ { return (**this).weight(); }
    auto& v() const _PHARE_ALL_FN_ { return (**this).v(); }

    Particle<dimension> copy() const _PHARE_ALL_FN_
    {
        return {(**this).weight(), (**this).charge(), (**this).iCell(), (**this).delta(),
                (**this).v()};
    }


    T* ts_particles_ptr;
    std::size_t idx = 0;
};




template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void TileSetSpan<Particles, impl>::sync() _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "TileSetSpan::sync()");

    auto view = *this;
    PHARE_WITH_MKN_GPU({
        mkn::gpu::GDLauncher{particles_.size()}([=] _PHARE_ALL_FN_() mutable {
            view.template sync_add_new<PHASE, type>();
            __syncthreads();
            view.template sync_rm_left<PHASE, type>();
        });
    })
}

template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type, typename... Args>
void TileSetSpan<Particles, impl>::sync(Args&&... args) _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "TileSetSpan::sync(stream)");

    auto view = *this;
    PHARE_WITH_MKN_GPU({
        auto const& [stream] = std::forward_as_tuple(args...);
        mkn::gpu::GDLauncher<true>{particles_.size()}.stream(stream, [=] _PHARE_ALL_FN_() mutable {
            view.template sync_add_new<PHASE, type>();
            __syncthreads();
            view.template sync_rm_left<PHASE, type>();
        });
    })
}


template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void TileSetSpan<Particles, impl>::sync_add_new() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU

    auto const& kidx = mkn::gpu::idx();
    auto& tile       = particles_[kidx];
    auto& real       = tile();
    auto const bix   = local_cell(tile.lower);

    auto const& n_gaps = gap_idx_(bix);
    {
        auto& gaps = gaps_(bix);
        thrust::sort(thrust::seq, gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
    }

    auto& left       = left_(bix);
    auto const& gaps = gaps_(bix);
    for (std::size_t i = 0; i < n_gaps; ++i)
    {
        auto const& gidx    = gaps[n_gaps - (1 + i)];
        auto const& newcell = local_tile_cell(real.iCell(gidx));
        auto const& cap     = cap_(newcell);
        auto& nparts        = (*particles_.at(newcell))();
        auto npidx          = nparts.size();

        while (true)
        {
            if (npidx >= cap)
                return;
            auto inc = npidx + 1;
            auto old = atomicCAS(nparts.size_address(), npidx, inc);
            if (npidx != old)
            {
                ++npidx;
                continue;
            }
            else
                break;
        }

        nparts.assign(real, gidx, npidx);
        ++left;
    }


#endif // PHARE_HAVE_MKN_GPU
}




template<typename Particles, std::uint8_t impl>
template<std::uint8_t PHASE, auto type>
void TileSetSpan<Particles, impl>::sync_rm_left() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU

    auto const& kidx = mkn::gpu::idx();
    auto& tile       = particles_[kidx];
    auto& real       = tile();
    auto const bix   = local_cell(tile.lower);
    auto const& gaps = gaps_(bix);
    auto& left       = left_(bix);
    auto& gaps_size  = gap_idx_(bix);
    add_into_(bix) -= left;
    while (left)
    {
        auto const& pidx = gaps[gaps_size - 1];
        real.assign(real.size() - 1, pidx);
        real.pop_back();
        --gaps_size;
        --left;
    }

#endif // PHARE_HAVE_MKN_GPU
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TILE_SET_HPP */
