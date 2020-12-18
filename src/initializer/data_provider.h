#ifndef DATA_PROVIDER_H
#define DATA_PROVIDER_H


#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include <functional>

#include "cppdict/include/dict.hpp"
#include "core/utilities/span.h"

namespace PHARE
{
namespace initializer
{
    template<typename Float, std::size_t dim>
    struct InitFunctionHelper
    {
    };

    template<>
    struct InitFunctionHelper<double, 1>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type)>;
    };

    template<>
    struct InitFunctionHelper<double, 2>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type, param_type)>;
    };

    template<>
    struct InitFunctionHelper<double, 3>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type, param_type, param_type)>;
    };



    template<>
    struct InitFunctionHelper<float, 1>
    {
        using return_type = std::shared_ptr<core::Span<float>>;
        using param_type  = std::vector<float> const&;
        using type        = std::function<return_type(param_type)>;
    };

    template<>
    struct InitFunctionHelper<float, 2>
    {
        using return_type = std::shared_ptr<core::Span<float>>;
        using param_type  = std::vector<float> const&;
        using type        = std::function<return_type(param_type, param_type)>;
    };

    template<>
    struct InitFunctionHelper<float, 3>
    {
        using return_type = std::shared_ptr<core::Span<float>>;
        using param_type  = std::vector<float> const&;
        using type        = std::function<return_type(param_type, param_type, param_type)>;
    };

    template<typename Float, std::size_t dim>
    using InitFunction = typename InitFunctionHelper<Float, dim>::type;



    using PHAREDict
        = cppdict::Dict<int, float, double, std::vector<double>, std::size_t,
                        std::optional<std::size_t>, std::string, InitFunction<float, 1>,
                        InitFunction<float, 2>, InitFunction<float, 3>, InitFunction<double, 1>,
                        InitFunction<double, 2>, InitFunction<double, 3>>;



    class PHAREDictHandler
    {
    public:
        static PHAREDictHandler& INSTANCE();

        void init() { phareDict = std::make_unique<PHAREDict>(); }

        void stop() { phareDict.release(); }

        auto& dict()
        {
            if (!phareDict)
                init();
            return *phareDict;
        }

    private:
        PHAREDictHandler() = default;
        std::unique_ptr<PHAREDict> phareDict;
    };


    class DataProvider
    {
    public:
        /**
         * @brief read will read the data from whatever source specific DataProvider classes
         * implement. readData() must be called before passing
         */
        virtual void read()     = 0;
        virtual ~DataProvider() = default;
    };

} // namespace initializer
} // namespace PHARE

#endif // DATA_PROVIDER_H
