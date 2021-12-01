#ifndef HYBRID_HYBRID_STATE_H
#define HYBRID_HYBRID_STATE_H


#include "core/models/physical_state.h"
#include "initializer/data_provider.h"

#include "core/utilities/algorithm.h"
#include "core/hybrid/hybrid_quantities.h"


#include <cstddef>
#include <sstream>
#include <string>
#include <utility>


namespace PHARE::core
{
template<typename HybridState_, typename GridLayout_>
struct HybridStateView
{
    static constexpr auto dimension = GridLayout_::dimension;

    using HybridState = HybridState_;
    using GridLayout  = GridLayout_;

    using ParticleArray_t = typename HybridState::ParticleArray_t;
    using VecField_t      = typename HybridState::VecField;
    using VecFieldView_t  = typename VecField_t::view_t;

    using Field_t     = typename VecField_t::field_type;
    using FieldView_t = typename Field_t::view_t;

    struct _EM_
    {
        VecFieldView_t E, B;
    };

    struct _IONS_
    {
        FieldView_t density;
        VecFieldView_t flux;
        ParticleArray_t* domain;
        ParticleArray_t* patch_ghost;
        ParticleArray_t* level_ghost;
        double const mass;
    };

    using Electromag_t = _EM_;
    using Ions_t       = _IONS_;


    HybridStateView(HybridState& state, GridLayout const& gridLayout)
        : layout{gridLayout}
        , electromag{state.electromag.E.view(), state.electromag.B.view()}
        , J{state.J.view()}
    {
        for (auto& pop : state.ions)
            ions.emplace_back(pop.density().view(), pop.flux().view(), &pop.domainParticles(),
                              &pop.patchGhostParticles(), &pop.levelGhostParticles(), pop.mass());
    }


    GridLayout layout;
    _EM_ electromag;
    std::vector<aggregate_adapter<_IONS_>> ions;
    VecFieldView_t J;
};

} // namespace PHARE::core

namespace PHARE
{
namespace core
{
    /**
     * @brief The HybridState class is a concrete implementation of a IPhysicalState.
     * It holds an Electromag, Ion and Electrons object manipulated by Hybrid concrete type of
     * ISolver
     */
    template<typename Electromag, typename Ions, typename Electrons>
    class HybridState : public IPhysicalState
    {
        using This = HybridState<Electromag, Ions, Electrons>;

    public:
        using VecField        = typename Electromag::vecfield_type;
        using ParticleArray_t = typename Ions::particle_array_type;

        template<typename GridLayout>
        using view_t = HybridStateView<This, GridLayout>;

        static constexpr auto dimension = Ions::dimension;

        HybridState(PHARE::initializer::PHAREDict const& dict)
            : electromag{dict["electromag"]}
            , ions{dict["ions"]}
            , J{"J", HybridQuantity::Vector::J}
            , electrons{dict["electrons"], ions, J}
        {
        }

        Electromag electromag;
        Ions ions;
        VecField J;
        Electrons electrons;

        std::string to_str()
        {
            std::stringstream ss;
            ss << "Hybrid State\n";
            ss << "------------------------------------\n";
            ss << core::to_str(ions);
            return ss.str();
        }


        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        bool isUsable() const { return electromag.isUsable() and ions.isUsable() && J.isUsable(); }



        bool isSettable() const
        {
            return electromag.isSettable() and ions.isSettable() && J.isSettable();
        }


        auto getCompileTimeResourcesUserList() const
        {
            return std::forward_as_tuple(electromag, ions, electrons);
        }

        auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(electromag, ions, electrons);
        }



        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------
    };



} // namespace core
} // namespace PHARE




#endif // PHARE_HYBRID_STATE_H
