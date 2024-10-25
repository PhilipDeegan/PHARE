#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/def.hpp"
#include "core/utilities/index/index.cpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"



namespace PHARE::core
{
enum class HyperMode : std::uint8_t { constant = 0, spatial, LAST };

template<typename GridLayout>
class Ohm_ref;

template<typename GridLayout>
class Ohm : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    explicit Ohm(PHARE::initializer::PHAREDict const& dict)
        : eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
        , hyper_mode{cppdict::get_value(dict, "hyper_mode", std::string{"constant"}) == "constant"
                         ? HyperMode::constant
                         : HyperMode::spatial}
    {
    }

    template<typename VecField, typename Field>
    void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                    VecField const& J, VecField& Enew) _PHARE_ALL_FN_
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ohm - GridLayout not set, cannot proceed to calculate ohm()");

        Ohm_ref<GridLayout>{*this->layout_, eta_, nu_, hyper_mode}(n, Ve, Pe, B, J, Enew);
    }


private:
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode;
};


template<typename GridLayout>
class Ohm_ref
{
    constexpr static auto dimension = GridLayout::dimension;
    // using This                      = Ohm_ref<GridLayout, VecField, Field>;

public:
    Ohm_ref(GridLayout const& layout_, double eta, double nu,
            HyperMode const hm = HyperMode::constant)
        : layout{layout_}
        , eta_{eta}
        , nu_{nu}
        , hyper_mode{hm}
    {
    }

    template<typename VecField, typename Field>
    void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                    VecField const& J, VecField& Enew) const _PHARE_ALL_FN_
    {
        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout.evalOnBox(
            Exnew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::X>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
        layout.evalOnBox(
            Eynew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::Y>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
        layout.evalOnBox(
            Eznew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::Z>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
    }

private:
    GridLayout layout; // GPU copy needs object not pointer!
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode;


    template<auto Tag, typename IJK, typename... Args>
    static void E_Eq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [n, Ve, Pe, B, J, E, self] = std::forward_as_tuple(args...);
        auto& Exyz                             = E(Tag);

        static_assert(Components::check<Tag>());

        Exyz(ijk) = self.template ideal_<Tag>(Ve, B, ijk)      //
                    + self.template pressure_<Tag>(n, Pe, ijk) //
                    + self.template resistive_<Tag>(J, ijk)    //
                    + self.template hyperresistive_<Tag>(J, B, n, ijk);
    }




    template<auto component, typename VecField>
    auto ideal1D_(VecField const& Ve, VecField const& B, MeshIndex<1> index) const _PHARE_ALL_FN_
    {
        if constexpr (component == Component::X)
        {
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);

            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }



    template<auto component, typename VecField>
    auto ideal2D_(VecField const& Ve, VecField const& B, MeshIndex<2> index) const _PHARE_ALL_FN_
    {
        if constexpr (component == Component::X)
        {
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);
            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }


        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }



    template<auto component, typename VecField>
    auto ideal3D_(VecField const& Ve, VecField const& B, MeshIndex<3> index) const _PHARE_ALL_FN_
    {
        if constexpr (component == Component::X)
        {
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);
            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }



    template<auto component, typename VecField>
    auto ideal_(VecField const& Ve, VecField const& B,
                MeshIndex<dimension> index) const _PHARE_ALL_FN_
    {
        if constexpr (dimension == 1)
            return ideal1D_<component>(Ve, B, index);
        if constexpr (dimension == 2)
            return ideal2D_<component>(Ve, B, index);
        if constexpr (dimension == 3)
            return ideal3D_<component>(Ve, B, index);
    }




    template<auto component, typename Field>
    auto pressure_(Field const& n, Field const& Pe, MeshIndex<dimension> index) const _PHARE_ALL_FN_
    {
        if constexpr (component == Component::X)
        {
            auto const nOnEx = GridLayout::project(n, index, GridLayout::momentsToEx());

            auto gradPOnEx = layout.template deriv<Direction::X>(Pe, index); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (component == Component::Y)
        {
            if constexpr (dimension >= 2)
            {
                auto const nOnEy = GridLayout::project(n, index, GridLayout::momentsToEy());

                auto gradPOnEy
                    = layout.template deriv<Direction::Y>(Pe, index); // TODO : issue 3391

                return -gradPOnEy / nOnEy;
            }
            else
            {
                return 0.;
            }
        }

        else if constexpr (component == Component::Z)
        {
            if constexpr (dimension >= 3)
            {
                auto const nOnEz = GridLayout::project(n, index, GridLayout::momentsToEz());

                auto gradPOnEz
                    = layout.template deriv<Direction::Z>(Pe, index); // TODO : issue 3391

                return -gradPOnEz / nOnEz;
            }
            else
            {
                return 0.;
            }
        }
    }




    template<auto component, typename VecField>
    auto resistive_(VecField const& J, MeshIndex<dimension> index) const _PHARE_ALL_FN_
    {
        auto const& Jxyx = J(component);

        if constexpr (component == Component::X)
        {
            auto const jxOnEx = GridLayout::project(Jxyx, index, GridLayout::JxToEx());
            return eta_ * jxOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const jyOnEy = GridLayout::project(Jxyx, index, GridLayout::JyToEy());
            return eta_ * jyOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const jzOnEz = GridLayout::project(Jxyx, index, GridLayout::JzToEz());
            return eta_ * jzOnEz;
        }
    }

    template<auto component, typename VecField, typename Field>
    auto hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                         MeshIndex<VecField::dimension> index) const _PHARE_ALL_FN_
    {
        // if compile error, fix this function
        static_assert(static_cast<std::underlying_type_t<HyperMode>>(HyperMode::LAST) == 2);

        if (hyper_mode == HyperMode::constant)
            return constant_hyperresistive_<component>(J, index);

        // else
        return spatial_hyperresistive_<component>(J, B, n, index);
    }


    template<auto component, typename VecField>
    auto constant_hyperresistive_(VecField const& J,
                                  MeshIndex<VecField::dimension> index) const _PHARE_ALL_FN_
    { // TODO : https://github.com/PHAREHUB/PHARE/issues/3
        return -nu_ * layout.laplacian(J(component), index);
    }


    template<auto component, typename VecField, typename Field>
    auto spatial_hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                 MeshIndex<VecField::dimension> index) const _PHARE_ALL_FN_
    { // TODO : https://github.com/PHAREHUB/PHARE/issues/3

        auto const lvlCoeff        = 1. / core::pow(4, layout.levelNumber());
        auto constexpr min_density = 0.1;

        auto computeHR = [&](auto BxProj, auto ByProj, auto BzProj, auto nProj) {
            auto const BxOnE = GridLayout::project(B(Component::X), index, BxProj);
            auto const ByOnE = GridLayout::project(B(Component::Y), index, ByProj);
            auto const BzOnE = GridLayout::project(B(Component::Z), index, BzProj);
            auto const nOnE  = GridLayout::project(n, index, nProj);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);

            return -nu_ * (b / (nOnE + min_density) + 1) * lvlCoeff
                   * layout.laplacian(J(component), index);
        };

        if constexpr (component == Component::X)
        {
            return computeHR(GridLayout::BxToEx(), GridLayout::ByToEx(), GridLayout::BzToEx(),
                             GridLayout::momentsToEx());
        }
        if constexpr (component == Component::Y)
        {
            return computeHR(GridLayout::BxToEy(), GridLayout::ByToEy(), GridLayout::BzToEy(),
                             GridLayout::momentsToEy());
        }
        if constexpr (component == Component::Z)
        {
            return computeHR(GridLayout::BxToEz(), GridLayout::ByToEz(), GridLayout::BzToEz(),
                             GridLayout::momentsToEz());
        }
    }
};


} // namespace PHARE::core
#endif
