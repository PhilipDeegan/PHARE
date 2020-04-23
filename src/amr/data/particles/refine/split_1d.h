#ifndef PHARE_SPLIT_1D_H
#define PHARE_SPLIT_1D_H
namespace PHARE::amr
{
template<>
struct SplitInnerSetter</*dim=*/1, /*nbrRefineParticles*/ 2>
{
    template<typename Weights, typename WeightsVal>
    static void set_weights(Weights& weights, WeightsVal weight_val)
    {
        for (auto& weight : weights)
            weight = weight_val;
    }

    template<typename Deltas, typename DeltaVal>
    static void set_deltas(Deltas& deltas, DeltaVal delta_val)
    {
        deltas[0][0] = -delta_val;
        deltas[1][0] = +delta_val;
    }
};

template<>
struct SplitInnerSetter</*dim=*/1, /*nbrRefineParticles*/ 3>
{
    template<typename Weights, typename WeightsVals>
    static void set_weights(Weights& weights, WeightsVals weight_vals)
    {
        weights[0] = weight_vals[0];
        weights[1] = weight_vals[1];
        weights[2] = weight_vals[1];
    }

    template<typename Deltas, typename DeltaVal>
    static void set_deltas(Deltas& deltas, DeltaVal delta_val)
    {
        deltas[0][0] = 0.0;
        deltas[1][0] = -delta_val;
        deltas[2][0] = +delta_val;
    }
};


/*************************************************************************************
  dim = 1
  interp = 1
  nbrRefineParticles = 2
*/
template<>
class Splitter<1, 1, 2> : public ASplitter<1, 1, 2>
{
public:
    Splitter()
        : ASplitter{weight_val_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_  = 0.551569;
    constexpr static float weight_val_ = 0.5;
};


/*************************************************************************
  dim = 1
  interp = 1
  nbrRefineParticles = 3
*/
template<>
class Splitter<1, 1, 3> : public ASplitter<1, 1, 3>
{
public:
    Splitter()
        : ASplitter{weight_vals_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_     = 1;
    constexpr static float weight_vals_[] = {0.5f, 0.25f};
};


/*************************************************************************
  dim = 1
  interp = 2
  nbrRefineParticles = 2
*/
template<>
class Splitter<1, 2, 2> : public ASplitter<1, 2, 2>
{
public:
    Splitter()
        : ASplitter{weight_val_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_  = 0.663959f;
    constexpr static float weight_val_ = 0.5;
};


/*************************************************************************
  dim = 1
  interp = 2
  nbrRefineParticles = 3
*/
template<>
class Splitter<1, 2, 3> : public ASplitter<1, 2, 3>
{
public:
    Splitter()
        : ASplitter{weight_vals_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_     = 1.112033f;
    constexpr static float weight_vals_[] = {0.468137, 0.265931};
};


/*************************************************************************
  dim = 1
  interp = 3
  nbrRefineParticles = 2
*/
template<>
class Splitter<1, 3, 2> : public ASplitter<1, 3, 2>
{
public:
    Splitter()
        : ASplitter{weight_val_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_  = 0.752399f;
    constexpr static float weight_val_ = 0.5;
};


/*************************************************************************
  dim = 1
  interp = 3
  nbrRefineParticles = 3
*/
template<>
class Splitter<1, 3, 3> : public ASplitter<1, 3, 3>
{
public:
    Splitter()
        : ASplitter{weight_vals_, delta_val_}
    {
    }

protected:
    constexpr static float delta_val_     = 1.275922;
    constexpr static float weight_vals_[] = {0.473943, 0.263028f};
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_1D_H*/