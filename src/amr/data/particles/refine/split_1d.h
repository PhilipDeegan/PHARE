#ifndef PHARE_SPLIT_1D_H
#define PHARE_SPLIT_1D_H
namespace PHARE::amr
{
template<size_t _interpOrder, size_t _refinedParticlesNbr>
class ASplitter_1d : public ASplitter</*dim=*/1, _refinedParticlesNbr>
{
protected:
    static constexpr size_t interpOrder = _interpOrder;
};

template<size_t _interpOrder, size_t _refinedParticlesNbr>
class Splitter_1d : public ASplitter_1d<_interpOrder, _refinedParticlesNbr>
{
    Splitter_1d() /* This class should never be instantiated */ = delete;
};

template<>
class Splitter_1d<1, 2> : public ASplitter_1d<1, 2>
{
public:
    using Super                                 = ASplitter_1d<1, 2>;
    static constexpr size_t dimension           = Super::dimension;
    static constexpr size_t interpOrder         = Super::interpOrder;
    static constexpr size_t refinedParticlesNbr = Super::refinedParticlesNbr;
    using Super::deltas_;
    using Super::weights_;

    Splitter_1d()
    {
        for (auto& weight : weights_)
            weight = weight_val_;

        deltas_[0][0] = -delta_val_;
        deltas_[1][0] = +delta_val_;
    }

protected:
    constexpr static float weight_val_ = 0.5;
    constexpr static float delta_val_  = 0.277f;
};


template<>
class Splitter_1d<2, 2> : public Splitter_1d<1, 2>
{
};


template<>
class Splitter_1d<3, 2> : public Splitter_1d<1, 2>
{
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_1D_H*/