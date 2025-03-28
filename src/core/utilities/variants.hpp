#ifndef PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP
#define PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP

#include "core/logger.hpp"
#include "core/utilities/types.hpp"

#include <vector>
#include <variant>

namespace PHARE::core
{
template<typename T>
auto ptr_or_null_for_type()
{
    using Ret = std::conditional_t<std::is_const_v<T>, void const*, void*>;
    return [](T& arg) mutable -> Ret { return &arg; };
}

template<typename... Ts>
struct overloads : Ts...
{
    using Ts::operator()...;
};

template<typename... Ts>
overloads(Ts&&...) -> overloads<std::decay_t<Ts>...>;


template<typename... Args>
auto constexpr _visit_ptr_overloads(std::tuple<Args...>*)
{
    return overloads{ptr_or_null_for_type<Args>()...,
                     [](auto&) mutable -> void* { return nullptr; },
                     [](auto const&) mutable -> void* { return nullptr; }};
}


template<typename T, typename... Ts>
struct unique : std::type_identity<T>
{
};

template<typename... Ts, typename U, typename... Us>
struct unique<std::tuple<Ts...>, U, Us...>
    : std::conditional_t<(std::is_same_v<U, Ts> || ...), unique<std::tuple<Ts...>, Us...>,
                         unique<std::tuple<Ts..., U>, Us...>>
{
};

template<typename... Ts>
using unique_tuple = typename unique<std::tuple<>, Ts...>::type;



template<typename... Args>
auto constexpr visit_ptr_overloads()
{
    return _visit_ptr_overloads((unique_tuple<Args...>*){});
}



template<typename Type, typename Variants>
auto& get_as_ref_or_throw(Variants& variants, std::size_t const start = 0)
{
    for (std::size_t idx = start; idx < variants.size(); ++idx)
        if (auto type = std::visit(visit_ptr_overloads<Type>(), variants[idx]))
            return *reinterpret_cast<Type*>(type);

    throw std::runtime_error("No element in variant for type");
}


// ARGS MUST BE IN THE SAME ORDER AS VARIANT LIST TYPES!!!!!
template<typename... Args, typename Variants>
auto get_as_tuple_or_throw(Variants& variants, std::size_t start = 0)
{
    using Tuple               = std::tuple<Args...>;
    auto constexpr tuple_size = std::tuple_size_v<Tuple>;

    auto ptr_or_null = visit_ptr_overloads<Args...>();

    auto pointer_tuple = for_N<tuple_size>([&](auto i) mutable {
        using Type = std::tuple_element_t<i, Tuple>;

        for (std::size_t idx = start; idx < variants.size(); ++idx)
            if (auto ptr = std::visit(ptr_or_null, variants[idx]))
            {
                ++start;
                return reinterpret_cast<Type*>(ptr);
            }
        return (Type*){nullptr};
    });

    for_N<tuple_size>([&](auto i) {
        if (std::get<i>(pointer_tuple) == nullptr)
            throw std::runtime_error("No element in variant for type");
    });

    return for_N<tuple_size, for_N_R_mode::forward_tuple>(
        [&](auto i) -> auto& { return *std::get<i>(pointer_tuple); });
}


} // namespace PHARE::core


#endif /*PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP*/
