#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE::core
{

/// @brief To indicate the side of a boundary location.
enum class Side { LOWER = -1, UPPER = 1 };

/*
 * Definitions for boundary types in 1d, 2d, and 3d: copy-pasted from SAMRAI's
 * appu/CartesianBoundaryDefines.h, with C-style enums turned into scoped enums.
 */

//@{
//! @name Definitions for boundary types in 1d, 2d, and 3d:
namespace Bdry
{
    enum class Type {
        UNDEFINED = -1,

        FACE3D = 1,
        EDGE3D = 2,
        NODE3D = 3,

        EDGE2D = 1,
        NODE2D = 2,

        NODE1D = 1
    };
}
//@}

/*
 * Definitions for boundary array sizes in 1d, 2d, or 3d:
 */

//@{
//! @name Definitions for boundary array sizes in 1d, 2d, or 3d:
int const NUM_1D_NODES = 2;

int const NUM_2D_EDGES = 4;
int const NUM_2D_NODES = 4;

int const NUM_3D_FACES = 6;
int const NUM_3D_EDGES = 12;
int const NUM_3D_NODES = 8;
//@}

/*
 * Definitions for Face, Edge, and Node boundary locations:
 *
 * Note that these definitions are used only for:
 * - Node boundary locations in 1d (XLO, XHI only), or
 * - Edge boundary locations in 2d (XLO, XHI, YLO, YHI only), or
 * - Face boundary locations in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary locations (see source code for more
//! information):
namespace BdryLoc
{
    enum class Type { XLO = 0, XHI = 1, YLO = 2, YHI = 3, ZLO = 4, ZHI = 5 };

    /// @brief Return the side of a boundary location.
    /// @param boundaryLoc The boundary location.
    /// @return The boundary side.
    constexpr Side side(Type boundaryLoc)
    {
        switch (boundaryLoc)
        {
            case Type::XLO:
            case Type::YLO:
            case Type::ZLO: return Side::LOWER; break;

            case Type::XHI:
            case Type::YHI:
            case Type::ZHI: return Side::UPPER; break;

            default: throw std::runtime_error("Invalid Type.");
        }
    };

    /// @brief Return the direction of a boundary location.
    /// @param boundaryLoc The boundary location.
    /// @return The boundary direction.
    constexpr Direction direction(Type boundaryLoc)
    {
        switch (boundaryLoc)
        {
            case Type::XLO:
            case Type::XHI: return Direction::X; break;

            case Type::YLO:
            case Type::YHI: return Direction::Y; break;

            case Type::ZLO:
            case Type::ZHI: return Direction::Z; break;

            default: throw std::runtime_error("Invalid Type.");
        }
    };
} // namespace BdryLoc
//@}

/*
 * Definitions for Node boundary locations in 2d:
 */

//@{
//! @name Definitions for Node boundary locations in 2d:
namespace NodeBdyLoc2D
{
    enum class Type { XLO_YLO = 0, XHI_YLO = 1, XLO_YHI = 2, XHI_YHI = 3 };
}
//@}

/*
 * Definitions for Edge boundary locations in 3d:
 */

//@{
//! @name Definitions for Edge boundary locations in 3d:
namespace EdgeBdyLoc3D
{
    enum class Type {
        XLO_YLO = 0,
        XHI_YLO = 1,
        XLO_YHI = 2,
        XHI_YHI = 3,
        XLO_ZLO = 4,
        XHI_ZLO = 5,
        XLO_ZHI = 6,
        XHI_ZHI = 7,
        YLO_ZLO = 8,
        YHI_ZLO = 9,
        YLO_ZHI = 10,
        YHI_ZHI = 11
    };
}
//@}

/*
 * Definitions for Node boundary locations in 3d:
 */

//@{
//! @name Definitions for Node boundary locations in 3d:
namespace NodeBdyLoc3D
{
    enum class Type {
        XLO_YLO_ZLO = 0,
        XHI_YLO_ZLO = 1,
        XLO_YHI_ZLO = 2,
        XHI_YHI_ZLO = 3,
        XLO_YLO_ZHI = 4,
        XHI_YLO_ZHI = 5,
        XLO_YHI_ZHI = 6,
        XHI_YHI_ZHI = 7
    };
}
//@}

/*
 * Definitions for Face, Edge, and Node boundary conditions:
 *
 * Note that FLOW, REFLECT, DIRICHLET and NEUMANN are used only for:
 * - Node boundary conditions in 1d, or
 * - Edge boundary conditions in 2d, or
 * - Face boundary conditions in 3d.
 *
 * Note that [X, Y, Z]FLOW, [X, Y, Z]REFLECT, [X, Y, Z]DIRICHLET, and
 * [X, Y, Z]NEUMANN are used only for:
 * - Node boundary conditions in 2d (X and Y cases only), or
 * - Edge and Node boundary conditions in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary conditions (see source code for more
//! information):
namespace BdryCond
{
    enum class Type {
        FLOW       = 90,
        REFLECT    = 91,
        DIRICHLET  = 92,
        NEUMANN    = 93,
        XFLOW      = 900,
        YFLOW      = 901,
        ZFLOW      = 902,
        XREFLECT   = 910,
        YREFLECT   = 911,
        ZREFLECT   = 912,
        XDIRICHLET = 920,
        YDIRICHLET = 921,
        ZDIRICHLET = 922,
        XNEUMANN   = 930,
        YNEUMANN   = 931,
        ZNEUMANN   = 932
    };
}
//@}


/**
 * @brief Returns the mirrored index of @p index with respect to a boundary.
 * @tparam side Whether we are reflecting across the Lower or Upper boundary.
 * @tparam centering The staggering of the data (Primal cells or Dual nodes).
 * @param index The directional index to be reflected.
 * @param boundaryLimitIndex The reference index representing the boundary location.
 * @return The reflected directional index.
 * @details The definition of `boundaryLimitIndex` depends on @p centering:
 * - **Primal:** The index of the first/last layer of cells *inside* the physical domain.
 * - **Dual:** The index of the mesh locations exactly *on* the boundary.
 */
template<Side side, QtyCentering centering>
NO_DISCARD inline constexpr std::uint32_t boundary_mirrored(std::uint32_t const index,
                                                            std::uint32_t const boundaryLimitIndex)
{
    int32_t constexpr s = static_cast<int32_t>(side);
    int32_t const i     = static_cast<int32_t>(index);
    int32_t const b     = static_cast<int32_t>(boundaryLimitIndex);

    if constexpr (centering == QtyCentering::primal)
    {
        return static_cast<std::uint32_t>(i - 2 * (i - b));
    }
    else // if constexpr (centering == QtyCentering::dual)
    {
        return static_cast<std::uint32_t>(i - 2 * (i - b) + s);
    };
}

/**
 * @brief Mirrors a multidimensional @p point across a boundary plane
 * @tparam dimension The number of spatial dimensions
 * @tparam direction The axis along which to mirror (X, Y, or Z)
 * @tparam side Upper or Lower boundary
 * @tparam centering Primal or Dual centering along @p direction
 * @param point The input point to be mirrored
 * @param boundaryLimitPoint Any point representing the boundary location, typically returned by
 * @c GridLayout ' @c physicalStartIndex or @c physicalEndIndex depending on the @p side value
 * @return A new Point with the mirrored coordinate in the @p direction axis
 */
template<std::size_t dimension, Direction direction, Side side, QtyCentering centering>
inline constexpr Point<std::uint32_t, dimension>
boundary_mirrored(Point<std::uint32_t, dimension> const point,
                  Point<std::uint32_t, dimension> const boundaryLimitPoint)
{
    Point<std::uint32_t, dimension> symmetricPoint = point;

    constexpr std::size_t iDir = static_cast<std::size_t>(direction);

    symmetricPoint[iDir]
        = boundary_mirrored<side, centering>(point[iDir], boundaryLimitPoint[iDir]);

    return symmetricPoint;
}

} // namespace PHARE::core

#endif /* PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP */