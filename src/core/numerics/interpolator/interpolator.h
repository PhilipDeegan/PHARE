#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H


#include <array>
#include <cstddef>
#include <functional>

#include "core/data/grid/gridlayout.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/point/point.h"
#include "core/utilities/range/ranges.h"

#include "core/def.h"
#include "core/logger.h"
#include "core/operators.h"

namespace PHARE::core
{
//! return the size of the index and weights arrays for
//! interpolation at a given order
// Number of points where the interpolator of order interpOrder
// will deposit mass/momentum. This is hence the size of the
// index and weight arrays for interpolation at a given order
constexpr int nbrPointsSupport(int interpOrder)
{
    return interpOrder + 1;
}



/** \brief the class Weight aims at computing the weight coefficient for
 *  interpolation at a specific order
 *
 *  This class assumes the interpolation order is known at compile-time
 *  thus there are three specialization for orders 1, 2 and 3.
 *
 *  the class has only one method called computeWeight that takes three arguments:
 *
 *  \param[in] normalized position is the particle position in the cell normalized by grid
 * spacing \param[in] startIndex first grid index where to interpolate the field \param[out]
 * weights contains the nbrPointsSupport weights calculated
 *
 *  these three parameters are given for a specific direction (x, y or z)
 */
template<std::size_t interpOrder>
class Weighter
{
};



/** \brief Specialization of Weight for first order interpolation
 */
template<>
class Weighter<1>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(1)>& weights)
    {
        weights[1] = normalizedPos - static_cast<double>(startIndex);
        weights[0] = 1. - weights[1];
    }

    static constexpr int interp_order = 1;
};


/** \brief specialization of Weighter for second order interpolation
 */
template<>
class Weighter<2>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(2)>& weights)
    {
        auto index = startIndex + 1;
        auto delta = static_cast<double>(index) - normalizedPos;
        double coef1, coef2, coef3;
        coef1 = 0.5 + delta;
        coef2 = delta;
        coef3 = 0.5 - delta;

        weights[0] = 0.5 * coef1 * coef1;
        weights[1] = 0.75 - coef2 * coef2;
        weights[2] = 0.5 * coef3 * coef3;
    }

    static constexpr int interp_order = 2;
};



/** \brief specialization of Weighter for third order interpolation
 */
template<>
class Weighter<3>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(3)>& weights)
    {
        constexpr double _4_over_3 = 4. / 3.;
        constexpr double _2_over_3 = 2. / 3.;

        auto index   = static_cast<double>(startIndex) - normalizedPos;
        double coef1 = 1. + 0.5 * index;
        double coef2 = index + 1;
        double coef3 = index + 2;
        double coef4 = 1. - 0.5 * (index + 3);

        double coef2_sq  = coef2 * coef2;
        double coef2_cub = coef2_sq * coef2;
        double coef3_sq  = coef3 * coef3;
        double coef3_cub = coef3_sq * coef3;

        weights[0] = _4_over_3 * coef1 * coef1 * coef1;
        weights[1] = _2_over_3 - coef2_sq - 0.5 * coef2_cub;
        weights[2] = _2_over_3 - coef3_sq + 0.5 * coef3_cub;
        weights[3] = _4_over_3 * coef4 * coef4 * coef4;
    }

    static constexpr int interp_order = 3;
};




//! Interpol performs the interpolation of a field using precomputed weights at
//! indices starting at startIndex. The class is templated by the Dimensionality
template<std::size_t dim>
class MeshToParticle
{
};


template<std::size_t dimdex, typename GridLayout, typename Quantity, Quantity quantity,
         typename IndexWeights>
auto static resolve_start_index_and_weights(IndexWeights const& indexWeights)
{
    auto constexpr centerings                              = GridLayout::centering(quantity);
    auto const& [d_starts, d_weights, p_starts, p_weights] = indexWeights;

    if constexpr (centerings[dimdex] == QtyCentering::primal)
        return std::forward_as_tuple(p_starts[dimdex], p_weights[dimdex]);
    else
        return std::forward_as_tuple(d_starts[dimdex], d_weights[dimdex]);
}

/** \brief specialization of Interpol for 1D interpolation
 */
template<>
class MeshToParticle<1>
{
public:
    /** Performs the 1D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field
     * \param[in] startIndex is the first of the nbrPointsSupport indices where to interpolate
     * the field \param[in] weights are the nbrPointsSupport weights used for the interpolation
     */
    template<typename GridLayout, typename Quantity, Quantity quantity, typename Field,
             typename IndexWeights>
    inline auto op(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = resolve_start_index_and_weights<0, GridLayout, Quantity, quantity>(indexWeights);

        auto const& order_size = xWeights.size();
        auto fieldAtParticle   = 0.;

        for (auto ik = 0u; ik < order_size; ++ik)
        {
            fieldAtParticle += field(xStartIndex + ik) * xWeights[ik];
        }
        return fieldAtParticle;
    }
};


/**\brief Specialization of Interpol for 2D interpolation
 */
template<>
class MeshToParticle<2>
{
public:
    /** Performs the 2D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each
     * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
     * interpolate the field in both directions \param[in] weights are the nbrPointsSupport
     * weights used for the interpolation in both directions
     */
    template<typename GridLayout, typename Quantity, Quantity quantity, typename Field,
             typename IndexWeights>
    inline auto op(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = resolve_start_index_and_weights<0, GridLayout, Quantity, quantity>(indexWeights);
        auto const& [yStartIndex, yWeights]
            = resolve_start_index_and_weights<1, GridLayout, Quantity, quantity>(indexWeights);

        auto const& order_size = xWeights.size();

        double fieldAtParticle = 0.;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            double Yinterp = 0.;
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                Yinterp += field(xStartIndex + ix, yStartIndex + iy) * yWeights[iy];
            }
            fieldAtParticle += Yinterp * xWeights[ix];
        }

        return fieldAtParticle;
    }
};



/** \brief Specialization of Interpol for 3D interpolation
 */
template<>
class MeshToParticle<3>
{
public:
    /** Performs the 3D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each
     * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
     * interpolate the field in the 3 directions \param[in] weights are the nbrPointsSupport
     * weights used for the interpolation in the 3 directions
     */
    template<typename GridLayout, typename Quantity, Quantity quantity, typename Field,
             typename IndexWeights>
    inline auto op(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = resolve_start_index_and_weights<0, GridLayout, Quantity, quantity>(indexWeights);
        auto const& [yStartIndex, yWeights]
            = resolve_start_index_and_weights<1, GridLayout, Quantity, quantity>(indexWeights);
        auto const& [zStartIndex, zWeights]
            = resolve_start_index_and_weights<2, GridLayout, Quantity, quantity>(indexWeights);

        auto const& order_size = xWeights.size();

        double fieldAtParticle = 0.;
        for (auto ix = 0u; ix < order_size; ++ix)
        {
            double Yinterp = 0.;
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                double Zinterp = 0.;
                for (auto iz = 0u; iz < order_size; ++iz)
                {
                    Zinterp += field(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)
                               * zWeights[iz];
                }
                Yinterp += Zinterp * yWeights[iy];
            }
            fieldAtParticle += Yinterp * xWeights[ix];
        }
        return fieldAtParticle;
    }
};




//! ParticleToMesh projects a particle density and flux to given grids
template<std::size_t dim, typename Operator = core::Operators<double>>
class ParticleToMesh
{
};



/** \brief specialization of ParticleToMesh for 1D interpolation
 */
template<typename Op>
class ParticleToMesh<1, Op>
{
public: /** Performs the 1D interpolation
         * \param[in] density is the field that will be interpolated from the particle Particle
         * \param[in] xFlux is the field that will be interpolated from the particle Particle
         * \param[in] yFlux is the field that will be interpolated from the particle Particle
         * \param[in] zFlux is the field that will be interpolated from the particle Particle
         * \param[in] fieldCentering is the centering (dual or primal) of the field in each
         * direction \param[in] particle is the single particle used for the interpolation of
         * density and flux \param[in] startIndex is the first index for which a particle will
         * contribute \param[in] weights is the arrays of weights for the associated index
         */
    template<typename Field, typename VecField, typename Particle, typename Indexes,
             typename Weights>
    inline void operator()(Field& density, VecField& flux, Particle const& particle,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux] = flux();
        auto const& [xStartIndex]         = startIndex;
        auto const& [xWeights]            = weights;
        auto const& order_size            = xWeights.size();

        auto const partRho   = particle.weight;
        auto const xPartFlux = particle.v[0] * particle.weight;
        auto const yPartFlux = particle.v[1] * particle.weight;
        auto const zPartFlux = particle.v[2] * particle.weight;

        for (auto ik = 0u; ik < order_size; ++ik)
        {
            Op{density(xStartIndex + ik)} += partRho * xWeights[ik] * coef;
            Op{xFlux(xStartIndex + ik)} += xPartFlux * xWeights[ik] * coef;
            Op{yFlux(xStartIndex + ik)} += yPartFlux * xWeights[ik] * coef;
            Op{zFlux(xStartIndex + ik)} += zPartFlux * xWeights[ik] * coef;
        }
    }
};




/** \brief specialization of ParticleToMesh for 2D interpolation
 */
template<typename Op>
class ParticleToMesh<2, Op>
{
public: /** Performs the 2D interpolation
         * \param[in] density is the field that will be interpolated from the particle Particle
         * \param[in] xFlux is the field that will be interpolated from the particle Particle
         * \param[in] yFlux is the field that will be interpolated from the particle Particle
         * \param[in] zFlux is the field that will be interpolated from the particle Particle
         * \param[in] fieldCentering is the centering (dual or primal) of the field in each
         * direction \param[in] particle is the single particle used for the interpolation of
         * density and flux \param[in] startIndex is the first index for which a particle will
         * contribute \param[in] weights is the arrays of weights for the associated index
         */
    template<typename Field, typename VecField, typename Particle, typename Indexes,
             typename Weights>
    inline void operator()(Field& density, VecField& flux, Particle const& particle,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]      = flux();
        auto const& [xStartIndex, yStartIndex] = startIndex;
        auto const& [xWeights, yWeights]       = weights;
        auto const& order_size                 = xWeights.size();

        auto const partRho   = particle.weight * coef;
        auto const xPartFlux = particle.v[0] * particle.weight * coef;
        auto const yPartFlux = particle.v[1] * particle.weight * coef;
        auto const zPartFlux = particle.v[2] * particle.weight * coef;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                auto x = xStartIndex + ix, y = yStartIndex + iy;

                Op{density(x, y)} += partRho * xWeights[ix] * yWeights[iy];
                Op{xFlux(x, y)} += xPartFlux * xWeights[ix] * yWeights[iy];
                Op{yFlux(x, y)} += yPartFlux * xWeights[ix] * yWeights[iy];
                Op{zFlux(x, y)} += zPartFlux * xWeights[ix] * yWeights[iy];
            }
        }
    }


    template<typename Particles, typename Field, typename VecField, typename Indexes,
             typename Weights>
    inline void op_0(Particles const& particles, std::size_t pi, Field& density, VecField& flux,
                     Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]      = flux();
        auto const& [xStartIndex, yStartIndex] = startIndex;
        auto const& [xWeights, yWeights]       = weights;
        auto const& order_size                 = xWeights.size();

        auto const& weight   = particles.weight[pi];
        auto const& v        = particles.v[pi];
        auto const partRho   = weight * coef;
        auto const xPartFlux = v[0] * weight * coef;
        auto const yPartFlux = v[1] * weight * coef;
        auto const zPartFlux = v[2] * weight * coef;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                auto x = xStartIndex + ix, y = yStartIndex + iy;

                Op{density(x, y)} += partRho * xWeights[ix] * yWeights[iy];
                Op{xFlux(x, y)} += xPartFlux * xWeights[ix] * yWeights[iy];
                Op{yFlux(x, y)} += yPartFlux * xWeights[ix] * yWeights[iy];
                Op{zFlux(x, y)} += zPartFlux * xWeights[ix] * yWeights[iy];
            }
        }
    }

    template<typename Particles, typename Field, typename VecField, typename IndexWeights>
    inline void op_1(Particles const& particles, std::size_t start, Field& density, VecField& flux,
                     IndexWeights const& indexWeights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]   = flux();
        auto const& [startIndexs, weightss] = indexWeights;

        for (std::size_t index = 0; index < weightss.size() and index < particles.size() - start;
             ++index)
        {
            auto const pi          = index + start;
            auto const& startIndex = startIndexs[index];
            auto const& weights    = weightss[index];

            auto const& [xStartIndex, yStartIndex] = startIndex;
            auto const& [xWeights, yWeights]       = weights;
            auto const& order_size                 = xWeights.size();

            auto const& weight   = particles.weight[pi];
            auto const& v        = particles.v[pi];
            auto const partRho   = weight * coef;
            auto const xPartFlux = v[0] * weight * coef;
            auto const yPartFlux = v[1] * weight * coef;
            auto const zPartFlux = v[2] * weight * coef;

            for (auto ix = 0u; ix < order_size; ++ix)
            {
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    auto x = xStartIndex + ix, y = yStartIndex + iy;

                    Op{density(x, y)} += partRho * xWeights[ix] * yWeights[iy];
                    Op{xFlux(x, y)} += xPartFlux * xWeights[ix] * yWeights[iy];
                    Op{yFlux(x, y)} += yPartFlux * xWeights[ix] * yWeights[iy];
                    Op{zFlux(x, y)} += zPartFlux * xWeights[ix] * yWeights[iy];
                }
            }
        }
    }


    template<typename Particles, typename Field, typename VecField, typename IndexWeights>
    inline void op_2(Particles const& particles, std::size_t start, Field& density, VecField& flux,
                     IndexWeights const& indexWeights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]   = flux();
        auto const& [startIndexs, weightss] = indexWeights;

        for (std::size_t index = 0; index < weightss.size() and index < particles.size() - start;
             ++index)
        {
            auto const pi          = index + start;
            auto const& startIndex = startIndexs[index];
            auto const& weights    = weightss[index];

            auto const& [xStartIndex, yStartIndex] = startIndex;
            auto const& [xWeights, yWeights]       = weights;
            auto const& order_size                 = xWeights.size();

            auto const& weight   = particles.weight[pi];
            auto const& v        = particles.v[pi];
            auto const partRho   = weight * coef;
            auto const xPartFlux = v[0] * weight * coef;
            auto const yPartFlux = v[1] * weight * coef;
            auto const zPartFlux = v[2] * weight * coef;

            for (auto ix = 0u; ix < order_size; ++ix)
                for (auto iy = 0u; iy < order_size; ++iy)
                    Op{density(xStartIndex + ix, yStartIndex + iy)}
                    += partRho * xWeights[ix] * yWeights[iy];

            for (auto ix = 0u; ix < order_size; ++ix)
                for (auto iy = 0u; iy < order_size; ++iy)
                    Op{xFlux(xStartIndex + ix, yStartIndex + iy)}
                    += xPartFlux * xWeights[ix] * yWeights[iy];

            for (auto ix = 0u; ix < order_size; ++ix)
                for (auto iy = 0u; iy < order_size; ++iy)
                    Op{yFlux(xStartIndex + ix, yStartIndex + iy)}
                    += yPartFlux * xWeights[ix] * yWeights[iy];

            for (auto ix = 0u; ix < order_size; ++ix)
                for (auto iy = 0u; iy < order_size; ++iy)
                    Op{zFlux(xStartIndex + ix, yStartIndex + iy)}
                    += zPartFlux * xWeights[ix] * yWeights[iy];
        }
    }


    // seems slightly faster than op_0 even if more ugly/coupled.
    template<typename Particles, typename Field, typename VecField, typename Indexes,
             typename Weights, typename IndexWeightsFn, typename Interpolator, typename GridLayout>
    inline void op_3(Particles const& particles, Field& density, VecField& flux,
                     Indexes const& startIndex, Weights const& weights, Interpolator& interpolator,
                     IndexWeightsFn fn, GridLayout const& layout, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]      = flux();
        auto const& [xStartIndex, yStartIndex] = startIndex;
        auto const& [xWeights, yWeights]       = weights;
        auto const& order_size                 = xWeights.size();

        for (std::size_t pi = 0; pi < particles.size(); ++pi)
        {
            fn(interpolator, layout, particles.iCell(pi), particles.delta[pi]);

            auto const& weight   = particles.weight[pi];
            auto const& v        = particles.v[pi];
            auto const partRho   = weight * coef;
            auto const xPartFlux = v[0] * weight * coef;
            auto const yPartFlux = v[1] * weight * coef;
            auto const zPartFlux = v[2] * weight * coef;

            for (auto ix = 0u; ix < order_size; ++ix)
            {
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    auto x = xStartIndex + ix, y = yStartIndex + iy;

                    Op{density(x, y)} += partRho * xWeights[ix] * yWeights[iy];
                    Op{xFlux(x, y)} += xPartFlux * xWeights[ix] * yWeights[iy];
                    Op{yFlux(x, y)} += yPartFlux * xWeights[ix] * yWeights[iy];
                    Op{zFlux(x, y)} += zPartFlux * xWeights[ix] * yWeights[iy];
                }
            }
        }
    }
};




/** \brief specialization of ParticleToMesh for 3D interpolation
 */
template<typename Op>
class ParticleToMesh<3, Op>
{
public: /** Performs the 3D interpolation
         * \param[in] density is the field that will be interpolated from the particle Particle
         * \param[in] xFlux is the field that will be interpolated from the particle Particle
         * \param[in] yFlux is the field that will be interpolated from the particle Particle
         * \param[in] zFlux is the field that will be interpolated from the particle Particle
         * \param[in] fieldCentering is the centering (dual or primal) of the field in each
         * direction \param[in] particle is the single particle used for the interpolation of
         * density and flux \param[in] startIndex is the first index for which a particle will
         * contribute \param[in] weights is the arrays of weights for the associated index
         */
    template<typename Field, typename VecField, typename Particle, typename Indexes,
             typename Weights>
    inline void operator()(Field& density, VecField& flux, Particle const& particle,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xFlux, yFlux, zFlux]                   = flux();
        auto const& [xStartIndex, yStartIndex, zStartIndex] = startIndex;
        auto const& [xWeights, yWeights, zWeights]          = weights;
        auto const& order_size                              = xWeights.size();

        auto const partRho   = particle.weight * coef;
        auto const xPartFlux = particle.v[0] * particle.weight * coef;
        auto const yPartFlux = particle.v[1] * particle.weight * coef;
        auto const zPartFlux = particle.v[2] * particle.weight * coef;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                for (auto iz = 0u; iz < order_size; ++iz)
                {
                    Op{density(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)}
                    += partRho * xWeights[ix] * yWeights[iy] * zWeights[iz];

                    Op{xFlux(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)}
                    += xPartFlux * xWeights[ix] * yWeights[iy] * zWeights[iz];

                    Op{yFlux(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)}
                    += yPartFlux * xWeights[ix] * yWeights[iy] * zWeights[iz];

                    Op{zFlux(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)}
                    += zPartFlux * xWeights[ix] * yWeights[iy] * zWeights[iz];
                }
            }
        }
    }
};




/** \brief Interpolator is used to perform particle-mesh interpolations using
 * 1st, 2nd or 3rd order interpolation in 1D, 2D or 3D, on a given layout.
 */
template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
class Interpolator : private Weighter<interpOrder>
{
    using This        = Interpolator<dim, interpOrder, atomic_ops>;
    using P2MOperator = core::Operators<double, atomic_ops>;

    using ParticleVecField = tuple_fixed_type<double, 3>;
    using ParticleEBs      = std::vector<tuple_fixed_type<ParticleVecField, 2>>;

    // this calculates the startIndex and the nbrPointsSupport() weights for
    // dual field interpolation and puts this at the corresponding location
    // in 'startIndex' and 'weights'. For dual fields, the normalizedPosition
    // is offseted compared to primal ones.
    template<typename CenteringT, CenteringT centering, typename GridLayout, typename ICell,
             typename Delta>
    auto indexAndWeights_(GridLayout const& layout, ICell const& iCell_, Delta const& delta)
    {
        // dual weights require -.5 to take the correct position weight
        auto constexpr dual_offset = .5;

        auto const& [startIndex_, weights_] = [&]() {
            if constexpr (centering == QtyCentering::dual)
                return std::forward_as_tuple(dual_startIndex_, dual_weights_);
            else
                return std::forward_as_tuple(primal_startIndex_, primal_weights_);
        }();

        auto iCell = layout.AMRToLocal(Point{iCell_});
        for (auto iDim = 0u; iDim < dimension; ++iDim)
        {
            startIndex_[iDim]
                = iCell[iDim] - computeStartLeftShift<CenteringT, centering>(delta[iDim]);

            double normalizedPos = iCell[iDim] + delta[iDim];

            if constexpr (centering == QtyCentering::dual)
                normalizedPos -= dual_offset;

            weightComputer_.computeWeight(normalizedPos, startIndex_[iDim], weights_[iDim]);
        }
    }

public:
    auto static constexpr interp_order = interpOrder;
    auto static constexpr dimension    = dim;
    /**\brief interpolate electromagnetic fields on all particles in the range
     *
     * For each particle :
     *  - The function first calculates the startIndex and weights for interpolation at
     * order InterpOrder and in dimension dim for dual and primal nodes
     *  - then it uses Interpol<> to calculate the interpolation of E and B components
     * onto the particle.
     */
    template<typename PartIterator, typename Electromag, typename GridLayout>
    inline void operator()(PartIterator begin, PartIterator end, Electromag const& Em,
                           GridLayout const& layout, ParticleEBs& particleEBs)
    {
        PHARE_LOG_SCOPE("Interpolator::operator()");

        using Scalar             = HybridQuantity::Scalar;
        auto const& [Ex, Ey, Ez] = Em.E();
        auto const& [Bx, By, Bz] = Em.B();

        // for each particle, first calculate the startIndex and weights for dual and
        // primal quantities. then, knowing the centering (primal or dual) of each
        // electromagnetic component, we use Interpol to actually perform the
        // interpolation. the trick here is that the StartIndex and weights have only been
        // calculated twice, and not for each E,B component.

        PHARE_LOG_START("MeshToParticle::operator()");
        std::size_t eb_idx = 0;
        for (auto currPart = begin; currPart != end; ++currPart)
        {
            auto& iCell = currPart->iCell;
            auto& delta = currPart->delta;
            indexAndWeights_<QtyCentering, QtyCentering::dual>(layout, iCell, delta);
            indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, iCell, delta);

            auto indexWeights = std::forward_as_tuple(dual_startIndex_, dual_weights_,
                                                      primal_startIndex_, primal_weights_);

            auto& [pE, pB]        = particleEBs[eb_idx++];
            auto& [pEx, pEy, pEz] = pE;
            auto& [pBx, pBy, pBz] = pB;

            pEx = meshToParticle_.template op<GridLayout, Scalar, Scalar::Ex>(Ex, indexWeights);
            pEy = meshToParticle_.template op<GridLayout, Scalar, Scalar::Ey>(Ey, indexWeights);
            pEz = meshToParticle_.template op<GridLayout, Scalar, Scalar::Ez>(Ez, indexWeights);
            pBx = meshToParticle_.template op<GridLayout, Scalar, Scalar::Bx>(Bx, indexWeights);
            pBy = meshToParticle_.template op<GridLayout, Scalar, Scalar::By>(By, indexWeights);
            pBz = meshToParticle_.template op<GridLayout, Scalar, Scalar::Bz>(Bz, indexWeights);
        }
        PHARE_LOG_STOP("MeshToParticle::operator()");
    }

    // container version of above
    template<typename Particles, typename Electromag, typename GridLayout>
    inline void meshToParticle(Particles& particles, Electromag const& Em, GridLayout const& layout,
                               ParticleEBs& particleEBs)
    {
        (*this)(particles.begin(), particles.end(), Em, layout, particleEBs);
    }



    /**\brief interpolate electromagnetic fields on all particles in the range
     *
     * For each particle :
     *  - The function first calculates the startIndex and weights for interpolation at
     * order InterpOrder and in dimension dim for dual and primal nodes
     *  - then it uses Interpol<> to calculate the interpolation of E and B components
     * onto the particle.
     */
    template<std::size_t version = 0, typename PartIterator, typename VecField, typename GridLayout,
             typename Field>
    inline void operator()(PartIterator begin, PartIterator end, Field& density, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
    {
        auto& startIndex_ = primal_startIndex_;
        auto& weights_    = primal_weights_;

        // for each particle, first calculate the startIndex and weights
        // for dual and primal quantities.
        // then, knowing the centering (primal or dual) of each electromagnetic
        // component, we use Interpol to actually perform the interpolation.
        // the trick here is that the StartIndex and weights have only been calculated
        // twice, and not for each E,B component.

        PHARE_LOG_START("ParticleToMesh::operator()");


        for (auto currPart = begin; currPart != end; ++currPart)
        {
            // TODO #3375
            indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, currPart->iCell,
                                                                 currPart->delta);

            particleToMesh_(density, flux, *currPart, startIndex_, weights_, coef);
        }
        PHARE_LOG_STOP("ParticleToMesh::operator()");
    }

    // container version of above
    template<std::size_t version = 0, typename ParticleRange, typename VecField,
             typename GridLayout, typename Field>
    inline void particleToMesh(ParticleRange const& particleRange, Field& density, VecField& flux,
                               GridLayout const& layout, double coef = 1.)
    {
        auto constexpr is_contiguous = ParticleRange::is_contiguous;

        if constexpr (!is_contiguous)
            this->operator()<version>(particleRange.begin(), particleRange.end(), density, flux,
                                      layout, coef);
        else
        {
            auto& particles   = particleRange.begin().particles;
            auto& startIndex_ = primal_startIndex_;
            auto& weights_    = primal_weights_;
            using Particles   = std::decay_t<decltype(particles)>;

            if constexpr (version == 0 || version == 3)
            {
                if constexpr (version == 0)
                    for (std::size_t pi = 0; pi < particles.size(); ++pi)
                    {
                        indexAndWeights_<QtyCentering, QtyCentering::primal>(
                            layout, particles.iCell(pi), particles.delta[pi]);

                        particleToMesh_.op_0(particles, pi, density, flux, startIndex_, weights_,
                                             coef);
                    }

                if constexpr (version == 3)
                    particleToMesh_.op_3(
                        particles, density, flux, startIndex_, weights_, *this,
                        std::mem_fn(
                            &This::template indexAndWeights_<QtyCentering, QtyCentering::primal,
                                                             GridLayout, std::array<int, dim>,
                                                             std::array<double, dim>>),
                        layout, coef);
            }
            else
            {
                auto constexpr batch = 1024 * 50;
                std::size_t offset   = 0;
                for (auto& range : ranges<Particles const, typename Particles::const_iterator>(
                         particles, batch, offset))
                {
                    auto indexWeightsArrayTuple
                        = get_indexAndWeights<QtyCentering, QtyCentering::primal, batch>(
                            layout, particles, offset);

                    if constexpr (version == 1)
                        particleToMesh_.op_1(particles, offset, density, flux,
                                             indexWeightsArrayTuple, coef);

                    if constexpr (version == 2)
                        particleToMesh_.op_2(particles, offset, density, flux,
                                             indexWeightsArrayTuple, coef);

                    offset += batch;
                }
            }
        }
    }


    /**
     * @brief Given a delta and an interpolation order, deduce which lower index to start
     * traversing from
     */
    template<typename CenteringT, CenteringT Centering>
    static int computeStartLeftShift([[maybe_unused]] double delta)
    {
        static_assert(interpOrder > 0 and interpOrder < 4);

        // If this is no longer true, it should be handled here via if constexpr/etc

        if constexpr (interpOrder == 1)
        {
            if constexpr (Centering == QtyCentering::primal)
                return 0;
            else
                return (delta < .5 ? 1 : 0);
        }

        else if constexpr (interpOrder == 2)
        {
            if constexpr (Centering == QtyCentering::primal)
                return (delta < .5 ? 1 : 0);
            else
                return 1;
        }

        else if constexpr (interpOrder == 3)
        {
            if constexpr (Centering == QtyCentering::primal)
                return 1;
            else
                return (delta < .5 ? 2 : 1);
        }
    }


private:
    static_assert(dimension <= 3 && dimension > 0 && interpOrder >= 1 && interpOrder <= 3, "error");

    template<typename T, std::size_t size>
    using Array = std::array<T, size>;

    // maybe could be std::uint8_t?
    using Starts  = Array<uint16_t, dimension>;
    using Weights = Array<Array<double, nbrPointsSupport(interpOrder)>, dimension>;

    Weighter<interpOrder> weightComputer_;
    MeshToParticle<dimension> meshToParticle_;
    ParticleToMesh<dimension, P2MOperator> particleToMesh_;

    Starts dual_startIndex_;
    Weights dual_weights_;

    Starts primal_startIndex_;
    Weights primal_weights_;


    template<typename CenteringT, CenteringT centering, std::size_t S = 50000, typename GridLayout,
             typename Particles>
    auto static get_indexAndWeights(GridLayout const& layout, Particles const& particles,
                                    std::size_t start)
    {
        // dual weights require -.5 to take the correct position weight
        auto constexpr dual_offset = .5;

        Weighter<interpOrder> weightComputer_;

        Array<Starts, S> startIndexs;
        Array<Weights, S> weights;

        auto left = particles.size() - start;
        auto end  = S > left ? left : S;
        for (std::size_t i = 0, pi = start; i < end; ++i, ++pi)
        {
            auto iCell        = layout.AMRToLocal(Point{particles.iCell(pi)});
            auto const& delta = particles.delta[pi];

            for (auto iDim = 0u; iDim < dimension; ++iDim)
            {
                startIndexs[i][iDim]
                    = iCell[iDim] - computeStartLeftShift<CenteringT, centering>(delta[iDim]);

                double normalizedPos = iCell[iDim] + delta[iDim];

                if constexpr (centering == QtyCentering::dual)
                    normalizedPos -= dual_offset;

                weightComputer_.computeWeight(normalizedPos, startIndexs[i][iDim],
                                              weights[i][iDim]);
            }
        }

        return std::make_tuple(startIndexs, weights);
    }
};



} // namespace PHARE::core

#endif
