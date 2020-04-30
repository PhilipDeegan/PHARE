#ifndef PHARE_CORE_FUNC_H
#define PHARE_CORE_FUNC_H

#include <functional>

namespace PHARE::core
{
template<typename ReturnType, std::size_t dim, bool smt>
struct ScalarFunctionHelper
{
};

template<bool _smt>
struct ScalarFunctionHelper<double, 1, _smt>
{
    static constexpr bool smt = _smt;
    using type                = std::function<double(double)>;
};

template<bool _smt>
struct ScalarFunctionHelper<double, 2, _smt>
{
    static constexpr bool smt = _smt;
    using type                = std::function<double(double, double)>;
};

template<bool _smt>
struct ScalarFunctionHelper<double, 3, _smt>
{
    static constexpr bool smt = _smt;
    using type                = std::function<double(double, double, double)>;
};

} // namespace PHARE::core

#endif /*PHARE_CORE_FUNC_H*/