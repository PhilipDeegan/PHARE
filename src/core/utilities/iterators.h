
#ifndef PHARE_CORE_UTILITIES_ITERATORS_H
#define PHARE_CORE_UTILITIES_ITERATORS_H

namespace PHARE::core
{
template<typename T, typename Vector, bool is_const = std::is_const_v<T>>
struct wrapped_iterator : public std::conditional_t<is_const, typename Vector::const_iterator,
                                                    typename Vector::iterator>
{
    // auto constexpr static is_const = std::is_const_v<T>;
    // using value_type               = typename Vector::value_type;
    using outer_type = T;
    using iterator_t
        = std::conditional_t<is_const, typename Vector::const_iterator, typename Vector::iterator>;

    // using difference_type   = typename iterator_t::difference_type;
    // using value_type        = typename iterator_t::value_type;
    // using pointer           = typename iterator_t::pointer;
    // using reference         = typename iterator_t::reference;
    // using iterator_category = typename iterator_t::iterator_category;

    wrapped_iterator operator+(std::size_t i)
    {
        wrapped_iterator copy = *this;
        static_cast<iterator_t&>(copy) += i;
        return copy;
    }

    // wrapped_iterator& operator=(wrapped_iterator const& src) = default;

    // auto operator-(wrapped_iterator const& that) const
    // {
    //     auto copy = *this;
    //     copy.it -= std::distance(that.it, copy.it); // copy.it - that.it;
    //     return copy;
    // }
    // auto& operator--()
    // {
    //     --it;
    //     return *this;
    // }
    // auto& operator++()
    // {
    //     ++it;
    //     return *this;
    // }

    // bool operator==(wrapped_iterator const& that) const
    // {
    //     return this->container == that.container and this->it == that.it;
    // }
    // bool operator!=(wrapped_iterator const& that) const { return !(*this == that); }

    // auto operator-> () { return &(*it); }
    // auto operator-> () const { return &(*it); }

    // auto& operator*() { return *it; }
    // auto& operator*() const { return *it; }

    // auto& operator()() const { return container; }
    auto& operator()() { return container; }

    // auto& super() const { return static_cast<iterator_t&>(*this); }
    auto& super() { return static_cast<iterator_t&>(*this); }

    // iterator_t it;

    T* container; // might be const
};

} // namespace PHARE::core


namespace std
{
// using namespace PHARE::core;

// template<typename T, typename V>
// auto distance(wrapped_iterator<T, V> const& a, wrapped_iterator<T, V> const& b)
// {
//     return std::distance(a.it, b.it);
// }

} // namespace std


#endif /* PHARE_CORE_UTILITIES_ITERATORS_H */