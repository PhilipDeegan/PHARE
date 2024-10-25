#ifndef PHARE_CORE_UTILITIES_TYPES_HPP_
#define PHARE_CORE_UTILITIES_TYPES_HPP_



#include "core/def.hpp"
#include "core/logger.hpp"


#include <cassert>
#include <array>
#include <iomanip>
#include <optional>
#include <string>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>


#include "cppdict/include/dict.hpp"
#include "menum/include/magic_enum/magic_enum_switch.hpp"
#include "menum/include/magic_enum/magic_enum_utility.hpp"

namespace PHARE
{
namespace core
{
    enum class Basis { Magnetic, Cartesian };




    template<typename T>
    NO_DISCARD std::vector<T> arange(T start, T stop, T step = 1)
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }

    template<typename T>
    NO_DISCARD T norm(std::array<T, 3> vec)
    {
        auto squarreSum = std::inner_product(std::begin(vec), std::end(vec), std::begin(vec), 0.);
        return std::sqrt(squarreSum);
    }



    enum class Edge { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax };


    template<typename T> // this is so we can use struct {} initialization with
                         // shared_ptrs/forwarding
    struct aggregate_adapter : public T
    {
        template<class... Args>
        aggregate_adapter(Args&&... args)
            : T{std::forward<Args>(args)...}
        {
        }
    };

    template<typename... Args> // this is so we can specialize
    struct type_list           // templates with only the outter most type
    {
        using Tuple = std::tuple<Args...>; //
    };


    template<typename T, typename A = void>
    struct is_std_vector : std::false_type
    {
    };

    template<typename T>
    struct is_std_vector<std::vector<T>> : std::true_type
    {
    };

    template<typename T, typename A>
    struct is_std_vector<std::vector<T, A>> : std::true_type
    {
    };

    template<typename T>
    inline constexpr auto is_std_vector_v = is_std_vector<T>::value;

    template<typename T, std::size_t size>
    struct is_std_array : std::false_type
    {
    };

    template<typename T, std::size_t size>
    struct is_std_array<std::array<T, size>, size> : std::true_type
    {
    };

    template<typename T, std::size_t size>
    inline constexpr auto is_std_array_v = is_std_array<T, size>::value;



    template<typename Tuple, typename Func>
    void apply(Tuple tuple, Func func)
    {
        std::apply([&](auto&... args) { (func(args), ...); }, tuple);
    }

    template<typename Type, std::size_t Size> // std::array::fill is only constexpr in C++20 ffs
    constexpr void fill(Type value, std::array<Type, Size>& array)
    {
        for (std::size_t i = 0; i < Size; i++)
            array[i] = value;
    }

    template<std::size_t Constant>
    class StrongIntegralConstant
    {
    public:
        constexpr decltype(auto) operator()() const { return constant(); }

    protected:
        static constexpr std::integral_constant<std::size_t, Constant> constant{};
    };

    template<std::size_t Constant>
    class DimConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class InterpConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class RefinedParticlesConst : public StrongIntegralConstant<Constant>
    {
    };

    template<std::size_t To_Size, std::size_t From_Size, typename Type>
    NO_DISCARD constexpr std::array<Type, To_Size>
    sized_array(std::array<Type, From_Size> const& from)
    {
        static_assert(To_Size <= From_Size, "invalid sized_array Size template, too large");

        if constexpr (To_Size == From_Size)
            return from;

        std::array<Type, To_Size> to{};

        for (std::size_t i = 0; i < to.size(); i++)
            to[i] = from[i];

        return to;
    }


    template<std::size_t To_Size, typename... Args>
    NO_DISCARD constexpr auto as_sized_array(Args&&... args)
    {
        auto arr = std::array{std::forward<decltype(args)>(args)...};

        return sized_array<To_Size>(arr);
    }

    template<typename Type, std::size_t size>
    NO_DISCARD constexpr std::array<Type, size> ConstArray(Type val = 0)
    {
        std::array<Type, size> arr{};
        for (uint8_t i = 0; i < size; i++)
            arr[i] = val;
        return arr;
    }

    template<std::size_t size, typename FN>
    constexpr auto ConstArrayFrom(FN fn)
    {
        std::array<decltype(fn()), size> arr{};
        for (uint8_t i = 0; i < size; i++)
            arr[i] = fn();
        return arr;
    }

    template<typename Type>
    NO_DISCARD std::vector<Type> displacementFrom(std::vector<Type> const& input)
    {
        std::vector<Type> displs(input.size());
        Type off = 0;
        for (Type i = 0; i < static_cast<Type>(input.size()); i++)
        {
            displs[i] = off;
            off += input[i];
        }
        return displs;
    }

    template<typename T>
    struct StackVar
    {
        using value_type = T;


        template<typename... Args>
        StackVar(Args&&... args)
            : var(std::forward<Args...>(args...))
        {
        }

        T var;
    };

    template<typename T>
    NO_DISCARD T from_string(std::string const& s)
    {
        T t;
        std::stringstream ss(s);
        ss >> t;
        return t;
    }

    template<typename T>
    NO_DISCARD std::string to_string_with_precision(T const& a_value, std::size_t const len)
    {
        std::ostringstream out;
        out.precision(len);
        out << std::fixed << a_value;
        auto str = out.str();
        // last digit may be non 0 because of rounding
        // and we know at that decimal it should be so we force it
        str.replace(str.end() - 1, str.end(), 1, '0');
        return out.str();
    }

    template<typename T>
    NO_DISCARD auto to_string_fixed_width(T const& value, std::size_t const& precision,
                                          std::size_t const& width, char const& fill = '0')
    {
        std::ostringstream out;
        out.width(width);
        out.precision(precision);
        out << std::setfill(fill) << std::fixed << value;
        return out.str();
    }

    template<typename T, std::size_t... Is>
    constexpr auto gft_helper(std::index_sequence<Is...> const&&)
        -> decltype(std::make_tuple((Is, std::declval<T>())...));

    template<typename T, std::size_t N>
    constexpr auto get_fixed_tuple() -> decltype(gft_helper<T>(std::make_index_sequence<N>{}));

    template<typename T, std::size_t N>
    using tuple_fixed_type = decltype(get_fixed_tuple<T, N>());



    NO_DISCARD inline std::optional<std::string> get_env(std::string const& key)
    {
        if (const char* val = std::getenv(key.c_str()))
            return std::string{val};
        return std::nullopt;
    }

    NO_DISCARD inline std::string get_env(std::string const& key, std::string const& _default)
    {
        if (auto e = get_env(key))
            return *e;
        return _default;
    }


    template<typename T>
    NO_DISCARD inline T get_env_as(std::string const& key, T const& t)
    {
        if (auto e = get_env(key))
            return from_string<T>(*e);
        return t;
    }

} // namespace core
} // namespace PHARE

namespace PHARE::core
{
template<typename Container, typename Multiplies = typename Container::value_type>
NO_DISCARD Multiplies product(Container const& container, Multiplies mul = 1) _PHARE_ALL_FN_
{
    // std accumulate doesn't exist on GPU
    for (auto const& v : container)
        mul *= v;
    return mul;
}

template<typename Container, typename Return = typename Container::value_type>
NO_DISCARD Return sum(Container const& container, Return r = 0) _PHARE_HST_FN_
{
    return std::accumulate(container.begin(), container.end(), r);
}

template<typename Container, typename F>
NO_DISCARD auto sum_from(Container const& container, F fn)
{
    using value_type  = typename Container::value_type;
    using return_type = std::decay_t<std::invoke_result_t<F const&, value_type const&>>;
    return_type sum   = 0;
    for (auto const& el : container)
        sum += fn(el);
    return sum;
}


template<typename Type>
auto& deref(Type&& type) _PHARE_ALL_FN_
{
    if constexpr (std::is_pointer_v<std::decay_t<Type>>)
        return *type;
    else
        return type;
}


template<typename F>
NO_DISCARD auto generate_from(F&& f, std::size_t from, std::size_t to)
{
    assert(from <= to);
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    std::vector<value_type> v;
    std::size_t count = to - from;
    if (count > 0)
        v.reserve(count);
    for (std::size_t i = from; i < to; ++i)
        v.emplace_back(f(i));
    return v;
}

template<typename F>
NO_DISCARD auto generate_from(F&& f, std::size_t count)
{
    return generate_from(std::forward<F>(f), 0, count);
}


template<typename F, typename Container>
NO_DISCARD auto generate_from(F&& f, Container& container)
{
    using T          = typename std::decay_t<Container>::value_type;
    using value_type = std::decay_t<std::invoke_result_t<F, T&>>;
    std::vector<value_type> vec;
    vec.reserve(container.size());
    for (auto& v : container)
        vec.emplace_back(f(v));
    return vec;
}


template<typename F, typename Container>
NO_DISCARD auto generate_from(F&& f, Container const& container)
{
    using T          = typename std::decay_t<Container>::value_type;
    using value_type = std::decay_t<std::invoke_result_t<F, T const&>>;
    std::vector<value_type> vec;
    vec.reserve(container.size());
    for (auto& v : container)
        vec.emplace_back(f(v));
    return vec;
}


template<std::size_t Idx, typename F, typename Type, std::size_t Size>
auto constexpr generate_array__(F& f, std::array<Type, Size> const& arr)
{
    return f(arr[Idx]);
}
template<typename Type, std::size_t Size, typename F, std::size_t... Is>
auto constexpr generate_array_(F& f, std::array<Type, Size> const& arr,
                               std::integer_sequence<std::size_t, Is...>)
{
    return std::array{generate_array__<Is>(f, arr)...};
}
template<typename F, typename Type, std::size_t Size>
auto constexpr generate_from(F&& f, std::array<Type, Size> const& arr)
{
    return generate_array_(f, arr, std::make_integer_sequence<std::size_t, Size>{});
}
// try to dedupe const/non-const above/below
template<std::size_t Idx, typename F, typename Type, std::size_t Size>
NO_DISCARD auto constexpr generate_array__(F& f, std::array<Type, Size>&& arr)
{
    return f(arr[Idx]);
}
template<typename Type, std::size_t Size, typename F, std::size_t... Is>
NO_DISCARD auto constexpr generate_array_(F& f, std::array<Type, Size>&& arr,
                                          std::integer_sequence<std::size_t, Is...>)
{
    return std::array{generate_array__<Is>(f, arr)...};
}
template<typename F, typename Type, std::size_t Size>
NO_DISCARD auto constexpr generate_from(F&& f, std::array<Type, Size>&& arr)
{
    return generate_array_(f, arr, std::make_integer_sequence<std::size_t, Size>{});
}



// calls operator bool() or copies bool
auto constexpr static to_bool = [](auto const& v) { return bool{v}; };


template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr all(Container const& container, Fn fn = to_bool)
{
    return std::all_of(container.begin(), container.end(), fn);
}

template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr any(Container const& container, Fn fn = to_bool)
{
    return std::any_of(container.begin(), container.end(), fn);
}

template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr none(Container const& container, Fn fn = to_bool)
{
    return std::none_of(container.begin(), container.end(), fn);
}

auto constexpr static accessor = [](auto const& v, auto const& i) -> auto& { return v[i]; };
template<typename Container, typename Fn = decltype(accessor)>
NO_DISCARD auto constexpr max_from(Container const& container, Fn fn = accessor)
{
    assert(container.size());
    auto t = fn(container, 0);
    for (std::size_t i = 1; i < container.size(); ++i)
        if (auto const& e = fn(container, i); e > t)
            t = e;
    return t;
}

template<typename Container, typename Fn = decltype(accessor)>
NO_DISCARD auto constexpr min_from(Container const& container, Fn fn = accessor)
{
    assert(container.size());
    auto t = fn(container, 0);
    for (std::size_t i = 1; i < container.size(); ++i)
        if (auto const& e = fn(container, i); e < t)
            t = e;
    return t;
}


template<typename SignedInt, typename UnsignedInt>
bool diff_sign_int_equals(SignedInt const& i0, UnsignedInt const& i1)
{
    static_assert(std::is_unsigned_v<UnsignedInt>);
    static_assert(std::is_signed_v<SignedInt>);
    static_assert(sizeof(UnsignedInt) >= sizeof(SignedInt), "Bad int comparison!");
    if (i0 < 0)
        return false;
    return static_cast<UnsignedInt>(i0) == i1;
}


template<typename Int0, typename Int1>
bool int_equals(Int0 const& i0, Int1 const& i1)
{
    if constexpr (std::is_same_v<Int0, Int1>)
        return i0 == i1;
    else
    {
        if constexpr (std::is_unsigned_v<Int0> and std::is_signed_v<Int1>)
            return diff_sign_int_equals(i1, i0);
        if constexpr (std::is_unsigned_v<Int1> and std::is_signed_v<Int0>)
            return diff_sign_int_equals(i0, i1);
    }
    // reaching here == compiler error
}




void inline abort_if(bool b)
{
    if (b)
        std::abort();
}




auto inline float_equals(float const& a, float const& b, float diff = 1e-6)
{
    return std::abs(a - b) < diff;
}

auto inline float_equals(double const& a, double const& b, double diff = 1e-12)
{
#ifndef NDEBUG
    auto ret  = std::abs(a - b);
    auto pred = ret < diff;
    if (!pred)
    {
        PHARE_LOG_LINE_STR(to_string_with_precision(a, 22));
        PHARE_LOG_LINE_STR(to_string_with_precision(b, 22));
        PHARE_LOG_LINE_STR(to_string_with_precision(ret, 22));
    }
    return pred;
#else
    return std::abs(a - b) < diff;
#endif
}


auto inline float_not_equals(double const& a, double const& b, double diff = 1e-12)
{
#ifndef NDEBUG
    auto ret  = std::abs(a - b);
    auto pred = ret > diff;
    if (!pred)
    {
        PHARE_LOG_LINE_STR(to_string_with_precision(a, 22));
        PHARE_LOG_LINE_STR(to_string_with_precision(b, 22));
        PHARE_LOG_LINE_STR(to_string_with_precision(ret, 22));
    }
    return pred;
#else
    return std::abs(a - b) > diff;
#endif
}

template<typename T = std::uint16_t>
struct Apply
{
    template<T i>
    auto constexpr operator()()
    {
        return std::integral_constant<T, i>{};
    }
};

template<typename Apply, std::uint16_t... Is>
constexpr auto apply_N(Apply& f, std::integer_sequence<std::uint16_t, Is...> const&)
{
    if constexpr (!std::is_same_v<decltype(f.template operator()<0>()), void>)
        return std::make_tuple(f.template operator()<Is>()...);
    (f.template operator()<Is>(), ...);
}

template<std::uint16_t N, typename Apply>
constexpr auto apply_N(Apply&& f)
{
    return apply_N(f, std::make_integer_sequence<std::uint16_t, N>{});
}

enum class for_N_R_mode {
    make_tuple = 0,
    make_array,
    forward_tuple,
};

template<std::uint16_t N, auto M = for_N_R_mode::make_tuple, typename Fn>
constexpr auto for_N(Fn& fn) _PHARE_ALL_FN_
{
    /*  // how to use
        for_N<2>([](auto ic) {
            constexpr auto i = ic();
            // ...
        });
    */

    static_assert(std::is_same_v<decltype(M), for_N_R_mode>);
    using return_type
        = std::decay_t<std::invoke_result_t<Fn, std::integral_constant<std::uint16_t, 0>>>;
    constexpr bool returns = !std::is_same_v<return_type, void>;

    if constexpr (returns)
    {
        return std::apply(
            [&](auto... ics) {
                if constexpr (M == for_N_R_mode::make_tuple)
                    return std::make_tuple(fn(ics)...);
                else if constexpr (M == for_N_R_mode::make_array)
                    return std::array{fn(ics)...};
                else if constexpr (M == for_N_R_mode::forward_tuple)
                    return std::forward_as_tuple(fn(ics)...);
                else
                    throw std::runtime_error("unknown return mode");
            },
            apply_N<N>(Apply{}));
    }
    else
        std::apply([&](auto... ics) { (fn(ics), ...); }, apply_N<N>(Apply{}));
}


template<std::uint16_t N, auto M = for_N_R_mode::make_tuple, typename Fn>
constexpr auto for_N(Fn&& fn)
{
    return for_N<N, M>(fn);
}

template<std::uint16_t N, typename Fn>
NO_DISCARD constexpr auto for_N_all(Fn&& fn)
{
    return all(for_N<N, for_N_R_mode::make_array>(fn));
}

template<std::uint16_t N, typename Fn>
NO_DISCARD constexpr auto for_N_any(Fn&& fn)
{
    return any(for_N<N, for_N_R_mode::make_array>(fn));
}

template<std::size_t S>
bool inline float_equals(std::array<double, S> const& a, std::array<double, S> const& b,
                         double diff = 1e-15)
{
    return for_N_all<S>([&](auto i) { return float_equals(a[i], b[i], diff); });
}


template<typename Tuple, std::size_t... Is>
constexpr auto named_pair_seq_(Tuple, std::index_sequence<Is...> const&&)
    -> decltype(std::make_tuple(
        (Is, std::declval<std::pair<std::string, std::tuple_element_t<Is, Tuple>>>())...));

template<typename... Args>
auto constexpr named_pair_seq()
    -> decltype(named_pair_seq_(std::tuple<Args...>{},
                                std::make_index_sequence<sizeof...(Args)>{}));

template<typename... Args>
using NamedTuple = decltype(named_pair_seq<Args...>());

template<typename... Pairs>
auto make_named_tuple(Pairs&&... pairs)
{
    return std::make_tuple(pairs...);
}


std::uint64_t inline now_in_microseconds()
{
    return std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::system_clock::now().time_since_epoch())
        .count();
}

template<typename T, std::size_t S>
auto static constexpr as_tuple(std::array<T, S> const& arr)
{
    return for_N<S, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return arr[i]; });
};

template<std::string_view const&... Strs>
struct join_string_views
{
    static constexpr auto impl() noexcept
    {
        constexpr std::size_t len = (Strs.size() + ... + 0);
        std::array<char, len + 1> arr{};
        auto append = [i = 0, &arr](auto const& s) mutable {
            for (auto c : s)
                arr[i++] = c;
        };
        (append(Strs), ...);
        arr[len] = 0;
        return arr;
    }
    static constexpr auto arr = impl();
    static constexpr std::string_view value{arr.data(), arr.size() - 1};
};
template<std::string_view const&... Strs>
static constexpr auto join_string_views_v = join_string_views<Strs...>::value;

template<typename T, T t>
struct to_string_view
{
    static constexpr auto size()
    {
        if constexpr (std::is_signed_v<T>)
        {
            std::size_t len = t > 0 ? 1 : 2;
            if (t < -9 or t > 9)
                for (auto n = t; n; len++, n /= 10) {}
            return len;
        }
        else
        {
            std::size_t len = 1;
            if (t > 9)
                for (auto n = t; n; len++, n /= 10) {}
            return len;
        }
    }

    static constexpr auto impl() noexcept
    {
        constexpr std::size_t len = size();
        static_assert(len == 1);

        std::array<char, len + 1> arr{};
        if (t != 0)
        {
            std::uint16_t i = len - 1;
            for (auto n = t; n; n /= 10)
                arr[i--] = "0123456789"[(t < 0 ? -1 : 1) * (n % 10)];
            if (t < 0)
                arr[i--] = '-';
        }
        else
        {
            arr[0] = '0';
        }
        arr[len] = 0;
        return arr;
    }

    static constexpr auto arr = impl();
    static constexpr std::string_view value{arr.data(), arr.size() - 1};
};

template<typename T, T t>
static constexpr auto to_string_view_v = to_string_view<T, t>::value;




template<typename T0, typename T1, std::size_t... Is>
bool _array_equals(T0 const& a, T1 const& b, std::index_sequence<Is...> const&&) _PHARE_ALL_FN_
{
    return (... && (a[Is] == b[Is]));
}

template<typename T, std::size_t S> // array == doesn't exist on GPU!
bool array_equals(std::array<T, S> const& a, std::array<T, S> const& b) _PHARE_ALL_FN_
{
    return _array_equals(a, b, std::make_index_sequence<S>{});
}


template<typename As, typename T, std::size_t... Is>
auto _array_minus(T const& a, T const& b, std::index_sequence<Is...> const&&) _PHARE_ALL_FN_
{
    std::array<As, sizeof...(Is)> arr;
    ((arr[Is] = a[Is] - b[Is]), ...);
    return arr;
}

template<typename As, typename T, std::size_t S>
auto array_minus(std::array<T, S> const& a, std::array<T, S> const& b) _PHARE_ALL_FN_
{
    return _array_minus<As>(a, b, std::make_index_sequence<S>{});
}

template<typename T>
auto constexpr pow(T const v, std::uint16_t const power) _PHARE_ALL_FN_
{
    T out = 1;
    for (std::uint16_t i = 0; i < power; ++i)
        out *= v;
    return out;
}


// template<typename A, typename T...>
auto constexpr any_in(auto const a, auto&&... ts) _PHARE_ALL_FN_
{
    return ((a == ts) || ...);
}




} // namespace PHARE::core


#endif // PHARE_CORE_UTILITIES_TYPES_HPP_
