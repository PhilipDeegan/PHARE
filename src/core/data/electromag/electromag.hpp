#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP

#include <string>
#include <tuple>

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "initializer/data_provider.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"



namespace PHARE::core::basic
{


template<typename VecFieldT>
class Electromag
{
public:
    static constexpr std::size_t dimension = VecFieldT::dimension;

    template<typename Em>
    auto static flat(Em& em)
    {
        auto& [Ex, Ey, Ez] = em.E();
        auto& [Bx, By, Bz] = em.B();
        return std::forward_as_tuple(Ex, Ey, Ez, Bx, By, Bz);
    }

    auto flat() { return flat(*this); }
    auto flat() const { return flat(*this); }

    VecFieldT E;
    VecFieldT B;
};


} // namespace PHARE::core::basic


namespace PHARE
{
namespace core
{
    template<typename VecFieldT>
    class Electromag : public basic::Electromag<VecFieldT>
    {
        using Super = basic::Electromag<VecFieldT>;

    public:
        using Super::B;
        using Super::E;
        using vecfield_type                    = VecFieldT;
        static constexpr std::size_t dimension = VecFieldT::dimension;

        explicit Electromag(std::string name)
            : Super{{name + "_E", HybridQuantity::Vector::E},
                    {name + "_B", HybridQuantity::Vector::B}}
            , Binit_{}
        {
        }

        explicit Electromag(initializer::PHAREDict const& dict)
            : Super{{dict["name"].template to<std::string>() + "_"
                         + dict["electric"]["name"].template to<std::string>(),
                     HybridQuantity::Vector::E},
                    {dict["name"].template to<std::string>() + "_"
                         + dict["magnetic"]["name"].template to<std::string>(),
                     HybridQuantity::Vector::B}}
            , Binit_{dict["magnetic"]["initializer"]}
        {
        }


        template<typename GridLayout>
        void initialize(GridLayout const& layout)
        {
            Binit_.initialize(B, layout);
        }


        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const { return E.isUsable() && B.isUsable(); }

        NO_DISCARD bool isSettable() const { return E.isSettable() && B.isSettable(); }

        NO_DISCARD auto getCompileTimeResourcesViewList() const
        {
            return std::forward_as_tuple(E, B);
        }

        NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(E, B); }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        void copyData(Electromag const& source)
        {
            E.copyData(source.E);
            B.copyData(source.B);
        }


        auto operator()() { return std::forward_as_tuple(E, B); }
        auto operator()() const { return std::forward_as_tuple(E, B); }

        template<typename V>
        auto as(auto&& a, auto&&... args)
        {
            return V{a(E, args...), a(B, args...)};
            // std::array components{&E, &B};
            // return V{for_N<2, for_N_R_mode::make_array>(
            //     [&](auto i) { return a(components[i], args...); })};
        }

        template<typename V>
        auto as(auto&& a, auto&&... args) const
        {
            return V{a(E, args...), a(B, args...)};
            // std::array const components{&E, &B};
            // return V{for_N<2, for_N_R_mode::make_array>(
            //     [&](auto i) { return a(components[i], args...); })};
        }


        Super& super() { return *this; }
        Super const& super() const { return *this; }
        auto& operator*() { return super(); }
        auto& operator*() const { return super(); }


    private:
        VecFieldInitializer<dimension> Binit_;
    };
} // namespace core
} // namespace PHARE
#endif
