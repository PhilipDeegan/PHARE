#ifndef PHARE_SPLIT_2D_H
#define PHARE_SPLIT_2D_H
namespace PHARE::amr
{
template<size_t _interpOrder, size_t _refinedParticlesNbr>
class ASplitter_2d : public ASplitter</*dim=*/2, _refinedParticlesNbr>
{
protected:
    static constexpr size_t interpOrder = _interpOrder;
};

template<size_t _interpOrder, size_t _refinedParticlesNbr>
class Splitter_2d : public ASplitter_2d<_interpOrder, _refinedParticlesNbr>
{
    Splitter_2d() /* This class should never be instantiated */ = delete;
};

template<>
class Splitter_2d<1, 9> : public ASplitter_2d<1, 9>
{
public:
    using Super                                 = ASplitter_2d<1, 9>;
    static constexpr size_t dimension           = Super::dimension;
    static constexpr size_t interpOrder         = Super::interpOrder;
    static constexpr size_t refinedParticlesNbr = Super::refinedParticlesNbr;
    using Super::deltas_;
    using Super::weights_;

    Splitter_2d()
    {
        weights_[4] = weight_vals_[0];
        for (size_t i : {0, 2, 6, 8})
            weights_[i] = weight_vals_[2];
        for (size_t i : {1, 3, 5, 7})
            weights_[i] = weight_vals_[1];

        constexpr static size_t root = 3; // 3*3 grid
        for (size_t y = 0; y < root; y++)
            for (size_t x = 0; x < root; x++)
                deltas_[x + (y * root)] = {delta_vals_[x], delta_vals_[y]};
    }

protected:
    constexpr static std::array<float, 3> weight_vals_ = {{0.25, 0.125, 0.01625}};
    constexpr static std::array<float, 3> delta_vals_  = {{-.25, 0, .25}};
};


template<>
class Splitter_2d<2, 9> : public Splitter_2d<1, 9>
{
};


template<>
class Splitter_2d<3, 9> : public Splitter_2d<1, 9>
{
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_2D_H*/