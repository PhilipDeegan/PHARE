#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_HPP

#include "diagnostic/detail/h5typewriter.hpp"

#include "core/data/vecfield/vecfield_component.hpp"
#include "diagnostic/diagnostic_model_view.hpp"

namespace PHARE::diagnostic::h5
{
/*
 * Possible outputs
 *
 * /t#/pl#/p#/electromag_(B, E)/(x,y,z)
 */
template<typename H5Writer>
class ElectromagDiagnosticWriter : public H5TypeWriter<H5Writer>
{
public:
    using Super = H5TypeWriter<H5Writer>;
    using Super::checkCreateFileFor_;
    using Super::fileData_;
    using Super::h5Writer_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Super::writeGhostsAttr_;
    using Attributes = typename Super::Attributes;
    using GridLayout = typename H5Writer::GridLayout;
    using FloatType  = typename H5Writer::FloatType;

    ElectromagDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}

    void createFiles(DiagnosticProperties& diagnostic) override;

    void getDataSetInfo(DiagnosticProperties& diagnostic, std::size_t iLevel,
                        std::string const& patchID, Attributes& patchAttributes) override;

    void initDataSets(DiagnosticProperties& diagnostic,
                      std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, std::size_t maxLevel) override;

    void writeAttributes(
        DiagnosticProperties&, Attributes&,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&,
        std::size_t maxLevel) override;

private:
    auto isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                      std::string var)
    {
        return diagnostic.quantity == tree + var;
    };
};


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    std::string tree = "/";
    checkCreateFileFor_(diagnostic, fileData_, tree, "B", "E");
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                          std::size_t iLevel,
                                                          std::string const& patchID,
                                                          Attributes& patchAttributes)
{
    auto& h5Writer         = this->h5Writer_;
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto& [id, type] : core::Components::componentMap())
        {
            // highfive doesn't accept uint32 which ndarray.shape() is
            auto const& array_shape = vecF.getComponent(type).shape();
            attr[name][id]          = std::vector<std::size_t>(array_shape.data(),
                                                               array_shape.data() + array_shape.size());
            auto ghosts = GridLayout::nDNbrGhosts(vecF.getComponent(type).physicalQuantity());
            attr[name][id + "_ghosts_x"] = static_cast<std::size_t>(ghosts[0]);
            if constexpr (GridLayout::dimension > 1)
                attr[name][id + "_ghosts_y"] = static_cast<std::size_t>(ghosts[1]);
            if constexpr (GridLayout::dimension > 2)
                attr[name][id + "_ghosts_z"] = static_cast<std::size_t>(ghosts[2]);
        }
    };

    auto per_qty = [&](auto& qty) {
        if (isActiveDiag(diagnostic, "/", qty.name()))
            infoVF(qty, qty.name(), patchAttributes[lvlPatchID]);
    };
    auto per_qtys = [&](auto&&... args) { ((per_qty(args), ...)); };


    if constexpr (is_hybrid_model<typename Super::Model_t>)
        per_qtys(h5Writer.modelView().getE(), h5Writer.modelView().getB());
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap())
        {
            auto vFPath = path + "/" + key + "_" + id;
            h5Writer.template createDataSet<FloatType>(
                h5file, vFPath,
                null ? std::vector<std::size_t>(GridLayout::dimension, 0)
                     : attr[key][id].template to<std::vector<std::size_t>>());

            this->writeGhostsAttr_(
                h5file, vFPath, null ? 0 : attr[key][id + "_ghosts_x"].template to<std::size_t>(),
                null);

            if constexpr (GridLayout::dimension > 1)
                this->writeGhostsAttr_(
                    h5file, vFPath,
                    null ? 0 : attr[key][id + "_ghosts_y"].template to<std::size_t>(), null);

            if constexpr (GridLayout::dimension > 2)
                this->writeGhostsAttr_(
                    h5file, vFPath,
                    null ? 0 : attr[key][id + "_ghosts_z"].template to<std::size_t>(), null);
        }
    };

    auto initPatch = [&](auto& level, auto& attr, std::string patchID = "") {
        bool null = patchID.empty();
        std::string path{h5Writer.getPatchPathAddTimestamp(level, patchID)};

        auto per_qty = [&](auto& qty) {
            if (isActiveDiag(diagnostic, "/", qty.name()))
                initVF(path, attr, qty.name(), null);
        };
        auto per_qtys = [&](auto&&... args) { ((per_qty(args), ...)); };

        if constexpr (is_hybrid_model<typename Super::Model_t>)
            per_qtys(h5Writer.modelView().getE(), h5Writer.modelView().getB());
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}



template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto per_qty = [&](auto& qty) {
        if (isActiveDiag(diagnostic, "/", qty.name()))
            h5Writer.writeTensorFieldAsDataset(*fileData_.at(diagnostic.quantity),
                                               h5Writer.patchPath() + "/" + qty.name(), qty);
    };
    auto per_qtys = [&](auto&&... args) { ((per_qty(args), ...)); };

    if constexpr (is_hybrid_model<typename Super::Model_t>)
        per_qtys(h5Writer.modelView().getE(), h5Writer.modelView().getB());
}



template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    writeAttributes_(diagnostic, *fileData_.at(diagnostic.quantity), fileAttributes,
                     patchAttributes, maxLevel);
}


} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
