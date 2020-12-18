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
    template<typename Float, std::size_t dimension>
    class VecFieldInitializer
    {
    public:
        using InitializerFunction = initializer::InitFunction<Float, dimension>;

        VecFieldInitializer() = default;

        VecFieldInitializer(initializer::PHAREDict dict)
            : x_{dict["x_component"].template to<InitializerFunction>()}
            , y_{dict["y_component"].template to<InitializerFunction>()}
            , z_{dict["z_component"].template to<InitializerFunction>()}
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
                                  InitializerFunction const& init)
        {
            auto indices      = layout.ghostStartToEndIndices(field, /*includeEnd=*/true);
            auto const coords = layout.template indexesToCoordVectors</*WithField=*/true>(
                indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                    return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
                });

            // keep grid data alive
            std::shared_ptr<Span<Float>> gridPtr
                = std::apply([&](auto&... args) { return init(args...); }, coords);
            Span<Float>& grid = *gridPtr;

            for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
                std::apply([&](auto&... args) { field(args...) = grid[cell_idx]; },
                           indices[cell_idx]);
        }



        InitializerFunction x_;
        InitializerFunction y_;
        InitializerFunction z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_H
