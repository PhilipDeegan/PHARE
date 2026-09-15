#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP

#include <cassert>
#include <iterator>
#include <type_traits>


namespace PHARE
{
namespace core
{
    template<typename...>
    using tryToInstanciate = void;


    struct dummy
    {
        using type              = int;
        static type const value = 0;
    };


    /** \brief Iterable is satisfied by any type that can be passed to std::begin/std::end,
     * e.g. a Box or a BoxContainer
     */
    template<typename IterableCandidate>
    concept Iterable = requires(IterableCandidate c) {
        std::begin(c);
        std::end(c);
    };


    template<typename IterableCandidate>
    using is_iterable = std::enable_if_t<Iterable<IterableCandidate>, dummy::type>;


    template<typename IterableCandidate>
    constexpr static bool is_iterable_v = Iterable<IterableCandidate>;


    // Basic function
    template<typename T>
    constexpr void allsame(T)
    {
    }

    // Recursive function
    template<typename T, typename T2, typename... Ts,
             typename = std::enable_if_t<std::is_same<T, T2>::value>>
    constexpr void allsame([[maybe_unused]] T arg, T2 arg2, Ts... args)
    {
        allsame(arg2, args...);
    }


} // namespace core

} // namespace PHARE

#endif
