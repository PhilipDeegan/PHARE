#ifndef DATA_PROVIDER_H
#define DATA_PROVIDER_H

#include "core/data/func.h"
#include "cppdict/include/dict.hpp"

#include <map>
#include <string>
#include <variant>

namespace PHARE
{
namespace initializer
{
    template<std::size_t dim>
    using ScalarFunction = typename core::ScalarFunctionHelper<double, dim, /*smt = */ false>::type;


    // template<std::size_t dim>
    using PHAREDict
        = cppdict::Dict<int, double, std::vector<double>, std::size_t, std::optional<std::size_t>,
                        std::string, ScalarFunction<1>, ScalarFunction<2>, ScalarFunction<3>>;



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
