#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core // 3d only for now
{

template<std::size_t dim>
struct BoxRow
{
    Point<std::uint32_t, dim> start_;
    std::uint32_t size_;

    auto& local_start() const { return start_; }
    auto& size() const { return size_; }
};

template<std::size_t dim>
struct BoxRows
{
    Box<std::uint32_t, dim> box;
    std::uint32_t k;
    std::uint32_t j = box.lower[1];
    std::uint32_t s = box.upper[2] - box.lower[2] + 1;

    BoxRow<dim> iter{{k, j, box.lower[2]}, s};

    BoxRows& operator++()
    {
        iter.start_ = Point{k, ++j, box.lower[2]};
        return *this;
    }

    auto& operator*() const { return iter; }
    bool operator==(BoxRows const& that) const { return j == that.j; }
    bool operator!=(BoxRows const& that) const { return j != that.j; }

    void next()
    {
        j = box.lower[1];
        ++k;
        iter.start_ = Point{k, j, box.lower[2]};
    }
};

template<std::size_t dim>
struct BoxSlab
{
    using BoxRows_t = BoxRows<dim>;

    Box<std::uint32_t, dim> box;
    bool _end = false;
    BoxRows_t br{box, _end ? box.upper[0] : box.lower[0]};

    BoxRows_t begin() const
    {
        assert(!_end);
        return br;
    }
    BoxRows_t end() const { return {box, box.upper[0], box.upper[1] + 1, 0}; }

    void next() { br.next(); }
};

template<std::size_t dim>
struct BoxSlabber
{
    using BoxSlab_t = BoxSlab<dim>;

    Box<std::uint32_t, dim> box;
    bool _end       = false;
    std::uint32_t k = _end ? box.upper[0] + 1 : box.lower[0];

    BoxSlab_t slab{box, _end};

    bool operator==(BoxSlabber const& that) const { return k == that.k; }
    bool operator!=(BoxSlabber const& that) const { return k != that.k; }

    BoxSlabber& operator++()
    {
        slab.next();
        ++k;
        return *this;
    }
    auto& operator*() { return slab; }
};

template<std::size_t dim>
struct BoxSpan
{
    using BoxSlabber_t = BoxSlabber<dim>;

    Box<std::uint32_t, dim> box;

    BoxSlabber_t b{box};
    BoxSlabber_t e{box, true};

    auto begin() { return b; }
    auto end() { return e; }
    auto begin() const { return b; }
    auto end() const { return e; }
};


template<std::size_t dim>
auto make_box_span(Box<std::uint32_t, dim> const box)
{
    return BoxSpan<dim>{box};
}


} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
