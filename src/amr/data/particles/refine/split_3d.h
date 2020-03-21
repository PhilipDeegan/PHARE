#ifndef PHARE_SPLIT_3D_H
#define PHARE_SPLIT_3D_H
namespace PHARE::amr
{
template<size_t _interpOrder, size_t _refinedParticlesNbr>
class ASplitter_3d : public ASplitter</*dim=*/3, _refinedParticlesNbr>
{
protected:
    static constexpr size_t interpOrder = _interpOrder;
};

template<size_t _interpOrder, size_t _refinedParticlesNbr>
class Splitter_3d : public ASplitter_3d<_interpOrder, _refinedParticlesNbr>
{
    Splitter_3d() /* This class should never be instantiated */ = delete;
};

template<>
class Splitter_3d<1, 16> : public ASplitter_3d<1, 16>
{
public:
    using Super                                 = ASplitter_3d<1, 16>;
    static constexpr size_t dimension           = Super::dimension;
    static constexpr size_t interpOrder         = Super::interpOrder;
    static constexpr size_t refinedParticlesNbr = Super::refinedParticlesNbr;
    using Super::deltas_;
    using Super::weights_;

    Splitter_3d() {}

protected:
};

template<>
class Splitter_3d<2, 16> : public Splitter_3d<1, 16>
{
};


template<>
class Splitter_3d<3, 16> : public Splitter_3d<1, 16>
{
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_3D_H*/