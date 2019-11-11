#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H

#include "detail/highfive.h"

namespace PHARE
{
namespace diagnostic
{
    namespace h5
    {
        template<typename HighFiveDiagnostic>
        class ElectromagDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
        {
        public:
            using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
            using Attributes = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
            ElectromagDiagnosticWriter(HighFiveDiagnostic& hi5)
                : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
            {
            }
            void write(DiagnosticDAO&) override;
            void compute(DiagnosticDAO&) override {}
            void getDataSetInfo(std::string const& patchID, Attributes& patchAttributes) override;
            void initDataSets(std::vector<std::string> const& patchIDs,
                              Attributes& patchAttributes) override;
        };

        template<typename HighFiveDiagnostic>
        void ElectromagDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
            std::vector<std::string> const& patchIDs, Attributes& patchAttributes)
        {
        }

        template<typename HighFiveDiagnostic>
        void
        ElectromagDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(std::string const& patchID,
                                                                       Attributes& patchAttributes)
        {
        }

        template<typename HighFiveDiagnostic>
        void ElectromagDiagnosticWriter<HighFiveDiagnostic>::write([
            [maybe_unused]] DiagnosticDAO& diagnostic)
        {
            auto& outer = this->hi5_;

            for (auto* vecField : outer.modelView().getElectromagFields())
            {
                outer.writeVecFieldAsDataset(outer.patchPath() + "/" + vecField->name(), *vecField);
            }
        }
    } // namespace h5
} // namespace diagnostic
} // namespace PHARE

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
