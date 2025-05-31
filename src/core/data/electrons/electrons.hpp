#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/utilities/variants.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"

#include "initializer/data_provider.hpp"

#include "core/utilities/index/index.hpp"



#include <memory>
#include <variant>

namespace PHARE::core
{


template<typename Ions>
class StandardHybridElectronFluxComputer
{
public:
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ve_{"StandardHybridElectronFluxComputer_Ve", HybridQuantity::Vector::V}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return core::isUsable(ions_, J_, Ve_); }

    NO_DISCARD bool isSettable() const { return core::isSettable(ions_, J_, Ve_); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return ions_.density();
    }

    NO_DISCARD Field& density()
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return ions_.density();
    }

    NO_DISCARD VecField& velocity()
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return velocity because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return Ve_;
    }

    void computeDensity() {}

    void computeBulkVelocity(GridLayout const& layout)
    {
        auto const& Jx  = J_(Component::X);
        auto const& Jy  = J_(Component::Y);
        auto const& Jz  = J_(Component::Z);
        auto const& Vix = ions_.velocity()(Component::X);
        auto const& Viy = ions_.velocity()(Component::Y);
        auto const& Viz = ions_.velocity()(Component::Z);
        auto const& Ni  = ions_.density(); // gives the particle density, hence the electron density

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components are defined on primal
        layout.evalOnBox(Ni, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::project(Jx, arr, GridLayout::JxToMoments());
            auto const JyOnVy = GridLayout::project(Jy, arr, GridLayout::JyToMoments());
            auto const JzOnVz = GridLayout::project(Jz, arr, GridLayout::JzToMoments());

            Vex(arr) = Vix(arr) - JxOnVx / Ni(arr);
            Vey(arr) = Viy(arr) - JyOnVy / Ni(arr);
            Vez(arr) = Viz(arr) - JzOnVz / Ni(arr);
        });
    }

    auto& getIons() const { return ions_; }

private:
    Ions ions_;
    VecField J_;
    VecField Ve_;
};



template<typename FluxComputer>
class ElectronPressureClosure
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    ElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
        : flux_{flux}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
    {
    }

    virtual ~ElectronPressureClosure() {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD virtual bool isUsable() const { return core::isUsable(Pe_, flux_); }

    NO_DISCARD virtual bool isSettable() const { return core::isSettable(Pe_, flux_); }

    struct PressureProperty // needed?
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;


    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(flux_, Pe_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(flux_, Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    void virtual initialize(GridLayout const& layout) = 0;

    NO_DISCARD Field& pressure()
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - ! isothermal closure pressure not usable");
        return Pe_;
    }

    NO_DISCARD Field const& pressure() const
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - !! isothermal closure pressure not usable");
        return Pe_;
    }

    void virtual computePressure(GridLayout const& /*layout*/) = 0;


    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }

    NO_DISCARD auto& B() const { return get_as_ref_or_throw<VecField const>(resources); }
    NO_DISCARD auto& B() { return get_as_ref_or_throw<VecField>(resources); }

    NO_DISCARD auto& Te() const
    {
        auto const& field = get_as_ref_or_throw<Field const>(resources);
        assert(field.name() == "Te");
        return field;
    }

protected:
    FluxComputer flux_;
    Field Pe_;

    using Resources = std::variant<VecField, Field>;
    std::vector<Resources> resources;
};



template<typename FluxComputer>
class IsothermalElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux)
        : Super{dict, flux}
        , T0_{dict["pressure_closure"]["Te"].template to<double>()}
    {
    }

    void initialize(GridLayout const& layout) override {}

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure is not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * T0_; });
    }

private:
    double const T0_ = 0;
};



template<typename FluxComputer>
class PolytropicElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;
    using Super::resources;

    PolytropicElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux)
        : Super{dict, flux}
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
        Field Te_{"Te", HybridQuantity::Scalar::P};
        resources.emplace_back(Te_);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& layout) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto const& V_ = this->flux_.velocity();

        Field Te_ = get_as_ref_or_throw<Field>(resources);

        std::transform(this->Pe_.begin(), this->Pe_.end(), N_.begin(), Te_.begin(),
                       [](auto p, auto n) { return p / n; });



        layout.evalOnBox(Te_, [&](auto&... ijk) mutable { P_Eq_(V_, Te_, ijk...); });


        // std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
        //                [this](auto n) { return n * 0.1; });  // TODO utiliser autre chose que
        //                transform et aller chercher Tnew
    }

private:
    double const gamma_ = 5. / 3.;
    initializer::InitFunction<dim> Pe_init_;


    template<typename Field, typename VecField>
    void P_Eq_(VecField const& V, Field const& T, auto&... ijk) const
    {
        // Tnew(ijk...) = T(ijk...);
        T(ijk...) = advection_(V, T, ijk...);
        //              + compression_(V, T, {ijk...});
    }

    template<typename Field, typename VecField>
    auto advection_(VecField const& Ve, Field const& Te, auto&... ijk) const
    {
        if constexpr (dim == 1)
            return advection1D_(Ve, Te, ijk...);
        // if constexpr (dim == 2)
        //     return advection2D_(Ve, Te, index);
        // if constexpr (dim == 3)
        //     return advection3D_(Ve, Te, index);
        return 0.;
    }

    template<typename Field, typename VecField>
    auto advection1D_(VecField const& Ve, Field const& Te, auto&... ijk) const
    {
        auto const& Vx = Ve(Component::X);
        // auto gradTonMoment_X = GridLayout::template deriv<Direction::X>(Te, ijk...);
        return Vx(ijk...); //*gradTonMoment_X;
                           // return 0.0*Te(ijk...);
    }
};



template<typename FluxComputer>
class CGLElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;
    using Super::resources;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;

    CGLElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux,
                               VecField const& B)
        : Super{dict, flux}
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
        resources.emplace_back(B);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - CGL closure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [](auto n) { return n * 0.1; });
    }

private:
    double const gamma_ = 5. / 3.;
    initializer::InitFunction<dim> Pe_init_;
};



template<typename FluxComputer>
std::unique_ptr<ElectronPressureClosure<FluxComputer>>
ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux,
                               typename FluxComputer::VecField const& B)
{
    if (dict["pressure_closure"]["name"].template to<std::string>() == "isothermal")
    {
        return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "polytropic")
    {
        return std::make_unique<PolytropicElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "CGL")
    {
        return std::make_unique<CGLElectronPressureClosure<FluxComputer>>(dict, flux, B);
    }
    return nullptr;
}



template<typename FluxComputer>
class ElectronMomentModel
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux,
                        VecField const& B)
        : dict_{dict}
        , fluxComput_{flux}
        , B_{B}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict, flux, B)}
    {
        assert(pressureClosure_);
    }

    ElectronMomentModel(ElectronMomentModel const& self)
        : dict_{self.dict_}
        , fluxComput_{self.fluxComput_}
        , B_{self.B_}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict_, fluxComput_, B_)}
    {
        *pressureClosure_ = *self.pressureClosure_;
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        // return fluxComput_.isUsable() and pressureClosure_->isUsable();
        return fluxComput_.isUsable() /*and B_.isUsable()*/ and pressureClosure_->isUsable();
    }

    // NO_DISCARD bool isSettable() const { return fluxComput_.isSettable(); }  // TODO
    // pressureClosure needs also to be settable ?
    NO_DISCARD bool isSettable() const
    {
        return fluxComput_.isSettable() /*and B_.isSettable()*/ and pressureClosure_->isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    void initialize(GridLayout const& layout) { pressureClosure_->initialize(layout); }

    NO_DISCARD Field const& density() const { return fluxComput_.density(); }
    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }
    NO_DISCARD Field const& pressure() const { return pressureClosure_->pressure(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }
    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }
    NO_DISCARD Field& pressure() { return pressureClosure_->pressure(); }

    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_->computePressure(layout); }

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    VecField B_; // NOT USED HERE!
    std::unique_ptr<ElectronPressureClosure<FluxComputer>> pressureClosure_;
};



template<typename FluxComputer>
class Electrons : public LayoutHolder<typename FluxComputer::GridLayout>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    Electrons(initializer::PHAREDict const& dict, FluxComputer flux, VecField const& B)
        : dict_{dict}
        , momentModel_{dict, flux, B}
    {
    }

    Electrons(Electrons const& that) = default;

    void initialize(GridLayout const& layout) { momentModel_.initialize(layout); }

    void update(GridLayout const& layout)
    {
        if (isUsable())
        {
            momentModel_.computeDensity();
            momentModel_.computeBulkVelocity(layout);
            momentModel_.computePressure(layout);
        }
        else
            throw std::runtime_error("Error - Electron  is not usable");
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return momentModel_.isUsable(); }

    NO_DISCARD bool isSettable() const { return momentModel_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(momentModel_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(momentModel_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD Field const& density() const { return momentModel_.density(); }
    NO_DISCARD VecField const& velocity() const { return momentModel_.velocity(); }
    NO_DISCARD Field const& pressure() const { return momentModel_.pressure(); }

    NO_DISCARD Field& density() { return momentModel_.density(); }
    NO_DISCARD VecField& velocity() { return momentModel_.velocity(); }
    NO_DISCARD Field& pressure() { return momentModel_.pressure(); }

private:
    initializer::PHAREDict dict_;
    ElectronMomentModel<FluxComputer> momentModel_;
};

} // namespace PHARE::core


#endif
