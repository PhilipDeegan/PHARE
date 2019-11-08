#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H



template<typename HighFiveDiagnostic>
class FluidDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
{
public:
    using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
    using Attributes = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
    FluidDiagnosticWriter(HighFiveDiagnostic& hi5)
        : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override{};
    void getDataSetInfo(std::string const& patchID, Attributes& patchAttributes) override;
    void initDataSets(std::vector<std::string> const& patchIDs,
                      Attributes& patchAttributes) override;
};


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    std::vector<std::string> const& patchIDs, Attributes& patchAttributes)
{
}

template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(std::string const& patchID,
                                                               Attributes& patchAttributes)
{
}

template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::write([[maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& outer = this->hi5_;
    auto& ions  = outer.modelView().getIons();
    std::string path(outer.patchPath() + "/ions/");

    for (auto& pop : outer.modelView().getIons())
    {
        std::string popPath(path + "pop/" + pop.name() + "/");
        auto& density = pop.density();
        outer.writeNewDataSet(popPath + "density", density.data(), density.size());
        outer.writeVecFieldAsDataset(popPath + "flux", pop.flux());
    }

    auto& density = ions.density();
    outer.writeNewDataSet(path + "density", density.data(), density.size());
    outer.writeVecFieldAsDataset(path + "bulkVelocity", ions.velocity());
}


#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
