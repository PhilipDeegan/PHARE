#ifndef VECFIELD_INITIALIZER_H
#define VECFIELD_INITIALIZER_H

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"
#include "initializer/data_provider.h"

#include <array>

namespace PHARE
{
namespace core
{
    template<std::size_t dimension>
    class VecFieldInitializer
    {
    public:
        VecFieldInitializer() = default;

        VecFieldInitializer(initializer::PHAREDict dict)
            : x_{dict["x_component"].template to<initializer::InitFunction<dimension>>()}
            , y_{dict["y_component"].template to<initializer::InitFunction<dimension>>()}
            , z_{dict["z_component"].template to<initializer::InitFunction<dimension>>()}
        {
        }


        template<typename VecField, typename GridLayout>
        void initialize(VecField& v, GridLayout const& layout)
        {
            static_assert(GridLayout::dimension == VecField::dimension,
                          "dimension mismatch between vecfield and gridlayout");

            initializeComponent_(v.getComponent(Component::X), layout, x_);
            initializeComponent_(v.getComponent(Component::Y), layout, y_);
            initializeComponent_(v.getComponent(Component::Z), layout, z_);
        }

    private:
        template<typename Field, typename GridLayout>
        void initializeComponent_(Field& field, GridLayout const& layout,
                                  initializer::InitFunction<1> const& init)
        {
            auto const [ix0, ix1] = layout.ghostStartToEnd(field, Direction::X);
            auto const [x]        = layout.template ghostGrids</*ByField=*/true>(
                field,
                [](auto& gridLayout, auto& field_, auto const&... args) {
                    return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
                },
                /*plus=*/1);
            assert(x.size() == (ix1 - ix0 + 1));

            std::shared_ptr<Span<double>> gridPtr = init(x); // keep grid data alive
            Span<double>& grid                    = *gridPtr;
            std::size_t cell                      = 0;
            for (std::uint32_t ix = ix0; ix <= ix1; ++ix)
                field(ix) = grid[cell++];
        }
        template<typename Field, typename GridLayout>
        void initializeComponent_(Field& field, GridLayout const& layout,
                                  initializer::InitFunction<2> const& init)
        {
            auto const [ix0, ix1] = layout.ghostStartToEnd(field, Direction::X);
            auto const [iy0, iy1] = layout.ghostStartToEnd(field, Direction::Y);
            auto const [x, y]     = layout.template ghostGrids</*ByField=*/true>(
                field,
                [](auto& gridLayout, auto& field_, auto const&... args) {
                    return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
                },
                /*plus=*/1);
            assert(x.size() == (ix1 - ix0 + 1) * (iy1 - iy0 + 1));

            std::shared_ptr<Span<double>> gridPtr = init(x, y);
            Span<double>& grid                    = *gridPtr;
            std::size_t cell                      = 0;
            for (std::uint32_t ix = ix0; ix <= ix1; ++ix)
                for (std::uint32_t iy = iy0; iy <= iy1; ++iy)
                    field(ix, iy) = grid[cell++];
        }
        template<typename Field, typename GridLayout>
        void initializeComponent_(Field& field, GridLayout const& layout,
                                  initializer::InitFunction<3> const& init)
        {
            auto const [ix0, ix1] = layout.ghostStartToEnd(field, Direction::X);
            auto const [iy0, iy1] = layout.ghostStartToEnd(field, Direction::Y);
            auto const [iz0, iz1] = layout.ghostStartToEnd(field, Direction::Z);
            auto const [x, y, z]  = layout.template ghostGrids</*ByField=*/true>(
                field,
                [](auto& gridLayout, auto& field_, auto const&... args) {
                    return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
                },
                /*plus=*/1);
            assert(x.size() == (ix1 - ix0 + 1) * (iy1 - iy0 + 1) * (iz1 - iz0 + 1));

            std::shared_ptr<Span<double>> gridPtr = init(x, y, z);
            Span<double>& grid                    = *gridPtr;
            std::size_t cell                      = 0;
            for (std::uint32_t ix = ix0; ix <= ix1; ++ix)
                for (std::uint32_t iy = iy0; iy <= iy1; ++iy)
                    for (std::uint32_t iz = iz0; iz <= iz1; ++iz)
                        field(ix, iy, iz) = grid[cell++];
        }


        initializer::InitFunction<dimension> x_;
        initializer::InitFunction<dimension> y_;
        initializer::InitFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_H
