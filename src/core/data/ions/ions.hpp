#ifndef PHARE_IONS_HPP
#define PHARE_IONS_HPP


#include "core/def.hpp"
#include "core/utilities/algorithm.hpp"
#include "initializer/data_provider.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
// #include "particle_initializers/particle_initializer_factory.hpp"



#include <array>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <functional>


namespace PHARE
{
namespace core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions
    {
    public:
        using value_type                = IonPopulation;
        using field_type                = IonPopulation::field_type;
        using vecfield_type             = IonPopulation::vecfield_type;
        using Float                     = field_type::type;
        using tensorfield_type          = IonPopulation::tensorfield_type;
        using particle_array_type       = IonPopulation::particle_array_type;
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;


        Ions(Ions const&) = default;
        Ions(Ions&&)      = default;


        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : rho_{densityName(), HybridQuantity::Scalar::rho}
            , massDensity_{massDensityName(), HybridQuantity::Scalar::rho}
            , bulkVelocity_{"bulkVel", HybridQuantity::Vector::V}
            , populations_{generate_from(
                  [&dict](auto ipop) { //
                      return IonPopulation{dict["pop" + std::to_string(ipop)]};
                  },
                  dict["nbrPopulations"].template to<std::size_t>())}
            , sameMasses_{allSameMass_()}
            , momentumTensor_{"momentumTensor", HybridQuantity::Tensor::M}
        {
        }


        NO_DISCARD auto nbrPopulations() const { return populations_.size(); }
        NO_DISCARD auto size() const { return nbrPopulations(); }


        NO_DISCARD field_type const& density() const { return rho_; }
        NO_DISCARD field_type& density() { return rho_; }

        NO_DISCARD field_type const& massDensity() const
        {
            return sameMasses_ ? rho_ : massDensity_;
        }
        NO_DISCARD field_type& massDensity() { return sameMasses_ ? rho_ : massDensity_; }


        NO_DISCARD vecfield_type const& velocity() const { return bulkVelocity_; }
        NO_DISCARD vecfield_type& velocity() { return bulkVelocity_; }

        NO_DISCARD std::string static densityName() { return "rho"; }
        NO_DISCARD std::string static massDensityName() { return "massDensity"; }

        tensorfield_type const& momentumTensor() const { return momentumTensor_; }
        tensorfield_type& momentumTensor() { return momentumTensor_; }

        void computeDensity()
        {
            rho_.zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.density();
                core::transform(rho_, popDensity, rho_, std::plus<Float>{});
            }
        }
        void computeMassDensity()
        {
            massDensity_.zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.density();
                core::transform(
                    massDensity_, popDensity, massDensity_,
                    [&pop](auto const& n, auto const& pop_n) { return n + pop_n * pop.mass(); });
            }
        }


        void computeBulkVelocity()
        {
            // the bulk velocity is sum(pop_mass * pop_flux) / sum(pop_mass * pop_density)
            // if all populations have the same mass, this is equivalent to sum(pop_flux) /
            // sum(pop_density) sum(pop_density) is rho_ and already known by the time we get here.
            // sum(pop_mass * pop_flux) is massDensity_ and is computed by computeMassDensity() if
            // needed
            if (!sameMasses_)
                computeMassDensity();
            auto const& density = (sameMasses_) ? rho_ : massDensity_;

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                // account for mass only if populations have different masses
                std::function<Float(Float, Float)> plus = std::plus<Float>{};
                std::function<Float(Float, Float)> plusMass
                    = [&pop](Float const& v, Float const& f) { return v + f * pop.mass(); };

                auto&& [fx, fy, fz] = pop.flux()();

                core::transform(vx, fx, vx, sameMasses_ ? plus : plusMass);
                core::transform(vy, fy, vy, sameMasses_ ? plus : plusMass);
                core::transform(vz, fz, vz, sameMasses_ ? plus : plusMass);
            }

            core::transform(vx, density, vx, std::divides<Float>{});
            core::transform(vy, density, vy, std::divides<Float>{});
            core::transform(vz, density, vz, std::divides<Float>{});
        }


        void computeFullMomentumTensor()
        {
            assert(momentumTensor_.isUsable());
            momentumTensor_.zero();
            auto& mom = momentumTensor_;

            for (auto& pop : populations_)
            {
                assert(pop.momentumTensor().isUsable());
                auto& p_mom = pop.momentumTensor();
                auto p_mij  = p_mom.begin();
                auto mij    = mom.begin();
                for (; p_mij != p_mom.end(); ++p_mij, ++mij)
                    core::transform(*mij, *p_mij, *mij, std::plus<typename field_type::type>{});
            }
        }


        NO_DISCARD auto begin() { return std::begin(populations_); }
        NO_DISCARD auto end() { return std::end(populations_); }
        NO_DISCARD auto begin() const { return std::begin(populations_); }
        NO_DISCARD auto end() const { return std::end(populations_); }


        // in the following isUsable and isSettable the massDensity_ is not checked
        // because it is for internal use only so no object will ever need to access it.
        NO_DISCARD bool isUsable() const
        {
            bool usable
                = rho_.isUsable() and bulkVelocity_.isUsable() and momentumTensor_.isUsable();

            // if all populations have the same mass, we don't need the massDensity_
            usable &= (sameMasses_) ? true : massDensity_.isUsable();

            for (auto const& pop : populations_)
            {
                usable = usable and pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const
        {
            bool settable
                = rho_.isSettable() and bulkVelocity_.isSettable() and momentumTensor_.isSettable();

            // if all populations have the same mass, we don't need the massDensity_
            settable &= (sameMasses_) ? true : massDensity_.isSettable();

            for (auto const& pop : populations_)
            {
                settable = settable and pop.isSettable();
            }
            return settable;
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::vector<IonPopulation>& getRunTimeResourcesViewList()
        {
            return populations_;
        }

        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(bulkVelocity_, momentumTensor_, rho_, massDensity_);
        }



        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------


        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << "Ions\n";
            ss << "------------------------------------\n";
            ss << "number of populations  : " << nbrPopulations() << "\n";
            for (auto& pop : populations_)
                ss << core::to_str(pop);
            return ss.str();
        }


        NO_DISCARD bool sameMasses() const { return sameMasses_; }

        auto& operator[](std::size_t const i) { return populations_[i]; }
        auto& operator[](std::size_t const i) const { return populations_[i]; }

    private:
        bool allSameMass_() const
        {
            return all(populations_, [this](auto const& pop) { // arbitrary small diff
                return float_equals(pop.mass(), populations_.front().mass(), /*abs_tol=*/1e-10);
            });
        }




        field_type rho_;
        field_type massDensity_;
        vecfield_type bulkVelocity_;
        std::vector<IonPopulation> populations_;

        // note this is set at construction when all populations are added
        // in the future if some populations are dynamically created during the simulation
        // this should be updated accordingly
        bool sameMasses_{false};
        tensorfield_type momentumTensor_;
    };
} // namespace core
} // namespace PHARE



#endif
