#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_PARTICLE_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_PARTICLE_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"


// not sure this makes sense

namespace PHARE::diagnostic::vtkh5
{

template<typename H5Writer>
class ParticlesDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    ParticlesDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}
};



template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();

    VTKFileWriter initializer{diagnostic, this};

    auto const write_quantity = [&](auto& layout, auto const&, auto const iLevel) {
        //
    };

    modelView.onLevels(
        [&](auto const& lvl) {
            auto const ilvl = lvl.getLevelNumber();
            initializer.initFileLevel(ilvl);

            // auto boxes = modelView.localLevelBoxes(ilvl);
            // for (auto* vecField : this->h5Writer_.modelView().getElectromagFields())
            //     if (diagnostic.quantity == "/" + vecField->name())
            //         initializer.template initTensorFieldFileLevel<1>(ilvl, boxes);

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H */
