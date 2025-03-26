#ifndef PHARE_TILE_SET_HPP
#define PHARE_TILE_SET_HPP


#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "tile_set_mapper.hpp"

#include <array>
#include <tuple>
#include <string>
#include <utility>


namespace PHARE::core
{


template<typename Tile>
class TileSetView
{
public:
    using Box_t                     = Box<int, Tile::dimension>;
    static auto constexpr dimension = Tile::dimension;

    TileSetView(Box_t const& box, std::array<std::size_t, dimension> const& tile_size,
                std::array<std::uint32_t, dimension> const& shape, Tile* tiles,
                std::size_t tile_nbr, Tile** cells,
                std::array<std::uint32_t, dimension> const& nbr_cells)
        : box_{box}
        , tile_size_{tile_size}
        , shape_{shape}
        , tiles_{tiles, tile_nbr}
        , cells_{cells, nbr_cells}
    {
    }

    TileSetView(TileSetView const&) = default;
    TileSetView(TileSetView&&)      = default;


    TileSetView& operator=(TileSetView const&) = default;
    TileSetView& operator=(TileSetView&&)      = default;


    NO_DISCARD auto overlaped_with(Box_t const& box) const
    {
        std::vector<std::pair<bool, Tile const*>> overlaped;
        for (auto const& tile : tiles_)
        {
            auto overlap = box * tile;
            if (overlap)
            {
                auto complete_overlap = (*overlap).size() == tile.size();
                overlaped.emplace_back(complete_overlap, &tile);
            }
        }
        return overlaped;
    }

    NO_DISCARD auto overlaped_with(Box_t const& box)
    {
        std::vector<std::pair<bool, Tile*>> overlaped;
        for (auto& tile : tiles_)
        {
            auto overlap = box * tile;
            if (overlap)
            {
                auto complete_overlap = (*overlap).size() == tile.size();
                overlaped.emplace_back(complete_overlap, &tile);
            }
        }
        return overlaped;
    }


    NO_DISCARD auto shape() const { return shape_; }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return tiles_.size(); }

    NO_DISCARD auto begin() { return tiles_.begin(); }
    NO_DISCARD auto begin() const { return tiles_.begin(); }

    NO_DISCARD auto end() { return tiles_.end(); }
    NO_DISCARD auto end() const { return tiles_.end(); }

    NO_DISCARD auto data() { return tiles_.data(); }
    NO_DISCARD auto data() const { return tiles_.data(); }

    NO_DISCARD auto& operator[](std::size_t const i) _PHARE_ALL_FN_ { return tiles_[i]; }
    NO_DISCARD auto& operator[](std::size_t const i) const _PHARE_ALL_FN_ { return tiles_[i]; }

    template<typename... Index>
    NO_DISCARD auto at(Index... indexes) _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_SS("");
        assert(cells_(indexes...));
        return cells_(indexes...);
    }
    template<typename... Index>
    NO_DISCARD auto at(Index... indexes) const _PHARE_ALL_FN_
    {
        // PHARE_LOG_LINE_SS("");
        assert(cells_(indexes...));
        return cells_(indexes...);
    }

    template<typename Tile_t>
    void reset(std::size_t i, Tile_t& tile)
    {
        tiles_[i].reset(tile);
    }

    template<typename TileSet_t>
    void reset(TileSet_t& tileset)
    {
        for (std::size_t i = 0; i < tileset.size(); ++i)
            reset(i, tileset[i]);
    }

private:
    Box_t const box_;
    std::array<std::size_t, dimension> tile_size_;
    std::array<std::uint32_t, dimension> shape_;
    Span<Tile> tiles_;
    NdArrayView<dimension, Tile*> cells_;
};



template<typename Tile, auto alloc_mode = AllocatorMode::CPU>
class TileSet
{
    bool static constexpr c_order = true;

    template<typename T>
    using nd_array_t = NdArrayVector<Tile::dimension, T, c_order, alloc_mode>;

    using Point_t = Point<std::uint32_t, Tile::dimension>;

public:
    template<typename T>
    using vector_t = typename PHARE::Vector<T, alloc_mode, 1>::vector_t;

    using value_type                = Tile;
    using Box_t                     = Box<int, Tile::dimension>;
    static auto constexpr dimension = Tile::dimension;

    template<typename... Args>
    TileSet(Box_t const& box, std::array<std::size_t, dimension> const& tile_size, Args&&... args)
        : box_{box}
        , tile_size_{tile_size}
        , shape_{[&]() {
            std::array<std::uint32_t, dimension> s;
            auto bs = box.shape();
            for (auto i = 0u; i < dimension; ++i)
                s[i] = (bs[i] + tile_size_[i] - 1) / tile_size_[i];
            return s;
        }()}
        , cells_{box.shape().template toArray<std::uint32_t>()}
    {
        tiles_.reserve(product(shape_));
        consistent_tile_size_();
        make_tiles_(args...);
        tag_cells_();
    }

    template<typename... Args>
    TileSet(Box_t const& box, std::size_t const& tile_size, Args&&... args)
        : TileSet{box, ConstArray<std::size_t, dimension>(tile_size), args...}
    {
    }

    TileSet(TileSet const& that)
        : box_{that.box_}
        , tile_size_{that.tile_size_}
        , shape_{that.shape_}
        , cells_{that.cells_}
        , tiles_{that.tiles_}
    {
        tag_cells_();
    }

    TileSet(TileSet&&) = default;

    TileSet& operator=(TileSet const& that)
    {
        box_       = that.box_;
        tile_size_ = that.tile_size_;
        shape_     = that.shape_;
        cells_     = that.cells_;
        tiles_     = that.tiles_;
        tag_cells_();
        return *this;
    }



    template<bool strict = false>
    NO_DISCARD auto overlaped_with(Box_t const& box)
    {
        if constexpr (strict)
            return strict_overlap_(box);
        else
            return lose_overlap_(box);
    }

    template<bool strict = false>
    NO_DISCARD auto overlaped_with(Box_t const& box) const
    {
        {
            if constexpr (strict)
                return strict_overlap_(box);
            else
                return lose_overlap_(box);
        }
    }

    NO_DISCARD auto inner_tiles() const
    {
        std::vector<Tile const*> border;
        Box_t inner_box{box_};
        inner_box.shrink(tile_size_);
        return overlaped_with<true>(inner_box);
    }

    NO_DISCARD auto inner_tiles()
    {
        std::vector<Tile*> border;
        Box_t inner_box{box_};
        inner_box.shrink(tile_size_);
        return overlaped_with<true>(inner_box);
    }

    // would be faster to have a way to directly select incomplete overlaps
    // rather than taking all and filtering incompletes
    NO_DISCARD auto border_tiles() const
    {
        for (auto const ts : tile_size_)
        {
            assert(ts > 1);
        }

        std::vector<Tile const*> border;
        Box_t inner_box{box_};
        // this box intersects the first tile around the perimeter
        // but completely contains the inner tiles
        // we want all incomplete overlaps
        auto tile_size_minus_one = tile_size_;
        for (auto& ts : tile_size_minus_one)
            --ts;
        inner_box.shrink(tile_size_minus_one);

        auto overlaped = overlaped_with(inner_box);
        for (auto const& [complete, tile] : overlaped)
        {
            if (!complete)
                border.push_back(tile);
        }
        return border;
    }

    NO_DISCARD auto border_tiles()
    {
        for (auto const ts : tile_size_)
        {
            assert(ts > 1);
        }

        std::vector<Tile const*> border;
        Box_t inner_box{box_};
        // this box intersects the first tile around the perimeter
        // but completely contains the inner tiles
        // we want all incomplete overlaps
        auto tile_size_minus_one = tile_size_;
        for (auto& ts : tile_size_minus_one)
            --ts;
        inner_box.shrink(tile_size_minus_one);

        auto overlaped = overlaped_with(inner_box);
        for (auto const& [complete, tile] : overlaped)
        {
            if (!complete)
                border.push_back(tile);
        }
        return border;
    }

    NO_DISCARD auto shape() const { return shape_; }
    NO_DISCARD auto size() const { return tiles_.size(); }

    NO_DISCARD auto begin() { return tiles_.begin(); }
    NO_DISCARD auto begin() const { return tiles_.begin(); }

    NO_DISCARD auto end() { return tiles_.end(); }
    NO_DISCARD auto end() const { return tiles_.end(); }

    NO_DISCARD auto data() { return tiles_.data(); }
    NO_DISCARD auto data() const { return tiles_.data(); }

    NO_DISCARD auto& operator[](std::size_t i) { return tiles_[i]; }
    NO_DISCARD auto const& operator[](std::size_t i) const { return tiles_[i]; }

    NO_DISCARD auto box() const { return box_; }
    NO_DISCARD auto tile_size() const { return tile_size_; }

    NO_DISCARD auto& operator()() { return tiles_; }
    NO_DISCARD auto& operator()() const { return tiles_; }


    NO_DISCARD auto& at(Point<int, dimension> const& amr_point)
    {
        return cells_((amr_point - box_.lower).as_unsigned());
    }
    NO_DISCARD auto& at(Point<int, dimension> const& amr_point) const
    {
        return cells_((amr_point - box_.lower).as_unsigned());
    }


    template<typename... Index>
    NO_DISCARD auto& at(Index... indexes)
    {
        return cells_(indexes...);
    }
    template<typename... Index>
    NO_DISCARD auto& at(Index... indexes) const
    {
        return cells_(indexes...);
    }


    template<typename View_t = Tile>
    auto make_view() // const ?
    {
        return TileSetView<View_t>{box_,          tile_size_,    shape_,        tiles_.data(),
                                   tiles_.size(), cells_.data(), cells_.shape()};
    }

    void static build_links(TileSet& tile_set)
    {
        auto const box_shape = tile_set.box_.shape().as_unsigned();

        if constexpr (dimension == 3)
        {
            // links[0] = {0, 0, 1};
            // links[1] = {0, 1, 0};
            // links[2] = {0, 1, 1};
            // links[3] = {1, 0, 0};
            // links[4] = {1, 0, 1};
            // links[5] = {1, 1, 0};
            // links[6] = {1, 1, 1};

            for (auto xi = 0u; xi < box_shape[0];)
            {
                Point const xip{xi, 0u, 0u};
                auto const xitile       = tile_set.cells_(xip);
                auto const xitile_shape = xitile->shape().as_unsigned();

                for (auto yi = 0u; yi < box_shape[1];)
                {
                    Point const yip{xi, yi, 0u};
                    auto const yitile       = tile_set.cells_(yip);
                    auto const yitile_shape = yitile->shape().as_unsigned();

                    for (auto zi = 0u; zi < box_shape[2];)
                    {
                        Point const zip{xi, yi, zi};
                        auto const zitile       = tile_set.cells_(zip);
                        auto const zitile_shape = zitile->shape().as_unsigned();
                        link_dim3(tile_set, box_shape, zip);
                        zi += zitile_shape[2];
                    }
                    yi += yitile_shape[1];
                }
                xi += xitile_shape[0];
            }
        }
    }

private:
    template<typename... Args>
    void static _link(Args&&... args)
    {
        auto const& [tile_set, tile, box_shape, point, idx] = std::forward_as_tuple(args...);
        if (point < box_shape)
            tile->link(idx) = tile_set.cells_(point);
    }

    template<typename... Args>
    void static link_dim3(Args&&... args)
    {
        auto const& [tile_set, box_shape, point] = std::forward_as_tuple(args...);
        auto tile                                = tile_set.cells_(point);
        auto const tile_shape                    = tile->shape().as_unsigned();

        {
            auto link0 = point;
            link0[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link0, 0);
        }
        {
            auto link1 = point;
            link1[1] += tile_shape[1];
            _link(tile_set, tile, box_shape, link1, 1);
        }
        {
            auto link2 = point;
            link2[1] += tile_shape[1];
            link2[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link2, 2);
        }
        {
            auto link3 = point;
            link3[0] += tile_shape[0];
            _link(tile_set, tile, box_shape, link3, 3);
        }
        {
            auto link4 = point;
            link4[0] += tile_shape[0];
            link4[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link4, 4);
        }
        {
            auto link5 = point;
            link5[0] += tile_shape[0];
            link5[1] += tile_shape[1];
            _link(tile_set, tile, box_shape, link5, 5);
        }

        _link(tile_set, tile, box_shape, point + tile_shape, 6);
    }




    template<typename TilePtr>
    static auto _lose_overlap(Box_t const& box, std::vector<Tile>& tiles)
    {
        std::vector<std::pair<bool, TilePtr>> overlaped;
        for (auto& tile : tiles)
        {
            auto overlap = box * tile;
            if (overlap)
            {
                auto complete_overlap = (*overlap).size() == tile.size();
                overlaped.emplace_back(complete_overlap, &tile);
            }
        }
        return overlaped;
    }
    template<typename TilePtr>
    NO_DISCARD static auto _strict_overlap(Box_t const& box, std::vector<Tile>& tiles)
    {
        std::vector<TilePtr> overlaped;
        for (auto& tile : tiles)
        {
            auto overlap = box * tile;
            if (overlap)
            {
                if (auto complete_overlap = (*overlap).size() == tile.size(); complete_overlap)
                {
                    overlaped.push_back(&tile);
                }
            }
        }
        return overlaped;
    }

    NO_DISCARD auto lose_overlap_(Box_t const& box) { return _lose_overlap<Tile*>(box, tiles_); }

    NO_DISCARD auto lose_overlap_(Box_t const& box) const
    {
        return _lose_overlap<Tile const*>(box, tiles_);
    }

    NO_DISCARD auto strict_overlap_(Box_t const& box)
    {
        return _strict_overlap<Tile*>(box, tiles_);
    }

    NO_DISCARD auto strict_overlap_(Box_t const& box) const
    {
        return _strict_overlap<Tile const*>(box, tiles_);
    }

    void consistent_tile_size_() const
    {
        auto const box_shape = box_.shape().template toArray<std::uint32_t>();
        for (auto idim = 0u; idim < dimension; ++idim)
        {
            if (box_shape[idim] < tile_size_[idim])
            {
                throw std::runtime_error("tile size larger than box size in dimension "
                                         + std::to_string(idim));
            }
        }
    }

    template<typename... Args>
    void make_tiles_(Args&&... args)
    {
        tile_set_make_tiles(*this, args...);
    }


    //! store the pointer to the tile associated with each cell
    void tag_cells_()
    {
        for (auto& tile : tiles_)
        {
            for (auto const& cell : tile)
            {
                // need to substract box lower to get
                // the local index of that cell in the NdArray
                cells_(cell - box_.lower) = &tile;
            }
        }
    }


    Box_t box_;
    std::array<std::size_t, dimension> tile_size_;
    std::array<std::uint32_t, dimension> shape_;
    nd_array_t<Tile*> cells_;
    vector_t<Tile> tiles_{};
};


template<typename F, typename Tile, auto alloc_mode = AllocatorMode::CPU>
auto& update_from(F f, TileSet<Tile, alloc_mode>& in)
{
    for (std::size_t i = 0; i < in.size(); ++i)
        in.data()[i] = f(i);
    return in;
}

template<auto alloc_mode = AllocatorMode::CPU, typename F, typename Tile, auto am1,
         typename... Args>
auto generate_from(F f, TileSet<Tile, am1>& in, Args&&... args)
{
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    TileSet<value_type, alloc_mode> ret{in.box(), in.tile_size(), args...};
    for (std::size_t i = 0; i < in.size(); ++i)
        ret.data()[i] = f(i);
    return ret;
}

template<auto alloc_mode = AllocatorMode::CPU, typename F, typename Tile, auto am1,
         typename... Args>
auto generate_from(F f, TileSet<Tile, am1> const& in, Args&&... args)
{
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    TileSet<value_type, alloc_mode> ret{in.box(), in.tile_size(), args...};
    for (std::size_t i = 0; i < in.size(); ++i)
        ret.data()[i] = f(i);
    return ret;
}

} // namespace PHARE::core


#endif
