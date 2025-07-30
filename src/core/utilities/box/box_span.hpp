#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core // 3d only for now
{

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct BoxRow
{
    using raw_value_type = Array_t::value_type;
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    std::uint32_t k;
    std::uint32_t j = box.lower[1];
    std::uint32_t s = box.upper[2] - box.lower[2] + 1;


    Span<value_type> row{&arr(k, j, box.lower[2]), s};

    BoxRow& operator++()
    {
        row = Span<value_type>{&arr(k, ++j, box.lower[2]), s};
        return *this;
    }
    auto& operator*() { return row; }

    bool operator==(BoxRow const& that) const { return j == that.j; }
    bool operator!=(BoxRow const& that) const { return j != that.j; }

    // double const* const last_domain_p_1 = &arr(box.upper);

    auto point() { return Point{k, j, box.lower[2]}; }

    void next()
    {
        j   = box.lower[1];
        row = Span<value_type>{&arr(++k, j, box.lower[2]), s};
    }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct BoxSlab
{
    using BoxRow_t = BoxRow<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end = false;
    BoxRow_t br{arr, box, _end ? box.upper[0] : box.lower[0]};

    BoxRow_t begin() const
    {
        assert(!_end);
        return br;
    }
    BoxRow_t end() const { return {arr, box, box.upper[0], box.upper[1] /*+ 1*/, 0, {0, 0}}; }

    void next() { br.next(); }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct BoxSlabber
{
    using BoxSlab_t = BoxSlab<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end       = false;
    std::uint32_t k = _end ? box.upper[0] /*+ 1 */ : box.lower[0];

    BoxSlab_t slab{arr, box, _end};

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

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct BoxSpan
{
    using BoxSlabber_t = BoxSlabber<Array_t>;

    Box<std::uint32_t, dim> box;
    Array_t& arr;

    BoxSlabber_t b{arr, box};
    BoxSlabber_t e{arr, box, true};

    auto begin() { return b; }
    auto end() { return e; }
    auto begin() const { return b; }
    auto end() const { return e; }
};


template<typename Array_t, std::size_t dim>
auto make_box_span(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return BoxSpan{box, arr};
}

template<typename Array_t, std::size_t dim>
auto make_const_box_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return BoxSpan<Array_t const>{box, arr};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
