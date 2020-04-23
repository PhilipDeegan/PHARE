#ifndef PHARE_SPLIT_2D_H
#define PHARE_SPLIT_2D_H
namespace PHARE::amr
{
template<>
struct SplitInnerSetter</*dim=*/2, /*nbrRefineParticles*/ 9>
{
    constexpr static size_t root = 3; // 3*3

    template<typename Weights, typename WeightsVal>
    static void set_weights(Weights& weights, WeightsVal& weight_vals)
    {
        weights[4] = weight_vals[0];
        for (size_t i : {0, 2, 6, 8})
            weights[i] = weight_vals[2];
        for (size_t i : {1, 3, 5, 7})
            weights[i] = weight_vals[1];
    }

    template<typename Deltas, typename DeltaVal>
    static void set_deltas(Deltas& deltas, DeltaVal& delta_vals)
    {
        for (size_t y = 0; y < root; y++)
            for (size_t yoff = y * root, x = 0; x < root; x++)
                deltas[x + yoff] = {delta_vals[x], delta_vals[y]};
    }
};


/*************************************************************************
  dim = 2
  interp = 1
  nbrRefineParticles = 9
*/
template<>
class Splitter<2, 1, 9> : public ASplitter<2, 1, 9>
{
public:
    Splitter()
        : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    constexpr static float delta_vals_[]  = {-.1, 0, .1};
    constexpr static float weight_vals_[] = {0.25, 0.125, 0.01625};
};


/*************************************************************************
  dim = 2
  interp = 2
  nbrRefineParticles = 9
*/
template<>
class Splitter<2, 2, 9> : public ASplitter<2, 2, 9>
{
public:
    Splitter()
    // : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    // constexpr static float delta_vals_[]  = {{-.1, 0, .1}};
    // constexpr static float weight_vals_[] = {{0.25, 0.125, 0.01625}};
};


/*************************************************************************
  dim = 2
  interp = 3
  nbrRefineParticles = 9
*/
template<>
class Splitter<2, 3, 9> : public ASplitter<2, 3, 9>
{
public:
    Splitter()
    // : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    // constexpr static float delta_vals_[]  = {{-.1, 0, .1}};
    // constexpr static float weight_vals_[] = {{0.25, 0.125, 0.01625}};
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_2D_H*/