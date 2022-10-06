#ifndef PHARE_FARADAY_AVX_HPP
#define PHARE_FARADAY_AVX_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "mkn/avx.hpp"

namespace PHARE::core
{
struct FaradayHelper
{
    static auto& I()
    {
        static FaradayHelper s;
        return s;
    }

    template<typename GridLayout>
    auto& set(GridLayout const& layout)
    {
        // Ez for primal in all directions - upper bound requirements
        auto curr_size = core::product(layout.allocSize(HybridQuantity::Scalar::Ez));
        if (derive_field.size() < curr_size)
            derive_field.resize(curr_size);
        return *this;
    }

    mkn::avx::Vector<double> derive_field;
};

template<typename GridLayout>
class FaradayVX : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");

        if (!(B.isUsable() && E.isUsable() && Bnew.isUsable()))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        this->dt_ = dt;

        BxEq_(B(Component::X), E, Bnew(Component::X));
        ByEq_(B(Component::Y), E, Bnew(Component::Y));
        BzEq_(B(Component::Z), E, Bnew(Component::Z));
    }


private:
    double dt_;


    template<typename VecField, typename Field, typename... Indexes>
    void BxEq_(Field const& Bx, VecField const& E, Field& Bxnew)
    {
        Bxnew.vector().assign(Bx.data(), Bx.data() + Bx.size());
        mkn::avx::Span<double> avx_BxNew{Bxnew};

        auto& avx_derive = FaradayHelper::I().derive_field;

        if constexpr (dimension == 2)
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), Bx.shape()};
            layout_->evalOnBox(Bxnew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::Y>(E(Component::Z), {ijk...});
            });
            avx_derive *= dt_;
            avx_BxNew -= avx_derive;
        }
        if constexpr (dimension == 3)
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), Bx.shape()};
            layout_->evalOnBox(Bxnew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::Y>(E(Component::Z), {ijk...});
            });
            avx_derive *= dt_;
            avx_BxNew -= avx_derive;
            layout_->evalOnBox(Bxnew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::Z>(E(Component::Y), {ijk...});
            });
            avx_derive *= dt_;
            avx_BxNew += avx_derive;
        }
    }

    template<typename VecField, typename Field, typename... Indexes>
    void ByEq_(Field const& By, VecField const& E, Field& Bynew)
    {
        Bynew.vector().assign(By.data(), By.data() + By.size());
        mkn::avx::Span<double> avx_ByNew{Bynew};

        auto& avx_derive = FaradayHelper::I().derive_field;

        if constexpr (dimension == 1 || dimension == 2)
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), By.shape()};
            layout_->evalOnBox(Bynew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::X>(E(Component::Z), {ijk...});
            });
            avx_derive *= dt_;
            avx_ByNew += avx_derive;
        }
        else
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), By.shape()};
            layout_->evalOnBox(Bynew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::Z>(E(Component::X), {ijk...});
            });
            avx_derive *= dt_;
            avx_ByNew -= avx_derive;
            layout_->evalOnBox(Bynew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::X>(E(Component::Z), {ijk...});
            });
            avx_derive *= dt_;
            avx_ByNew += avx_derive;
        }
    }

    template<typename VecField, typename Field, typename... Indexes>
    void BzEq_(Field const& Bz, VecField const& E, Field& Bznew)
    {
        Bznew.vector().assign(Bz.data(), Bz.data() + Bz.size());
        auto& avx_derive = FaradayHelper::I().derive_field;

        mkn::avx::Span<double> avx_Bznew{Bznew};

        if constexpr (dimension == 1)
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), Bz.shape()};
            layout_->evalOnBox(Bznew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::X>(E(Component::Y), {ijk...});
            });
            avx_derive *= dt_;
            avx_Bznew -= avx_derive;
        }

        else
        {
            NdArrayView<dimension, double> derive{avx_derive.data(), Bz.shape()};
            layout_->evalOnBox(Bznew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::X>(E(Component::Y), {ijk...});
            });
            avx_derive *= dt_;
            avx_Bznew -= avx_derive;
            layout_->evalOnBox(Bznew, [&](auto const&... ijk) {
                derive(ijk...) = layout_->template deriv<Direction::Y>(E(Component::X), {ijk...});
            });
            avx_derive *= dt_;
            avx_Bznew += avx_derive;
        }
    }
};

} // namespace PHARE::core


#endif
