#ifndef PHARE_SPLIT_3D_H
#define PHARE_SPLIT_3D_H
namespace PHARE::amr
{
template<>
struct SplitInnerSetter</*dim=*/3, /*nbrRefineParticles*/ 27>
{
    constexpr static size_t root = 3; // 3*3*3

    template<typename Weights, typename WeightsVal>
    static void set_weights(Weights& weights, WeightsVal& weight_vals)
    {
        // TODO
    }

    template<typename Deltas, typename DeltaVal>
    static void set_deltas(Deltas& deltas, DeltaVal& delta_vals)
    {
        for (size_t z = 0; z < root; z++)
            for (size_t zoff = z * root * root, y = 0; y < root; y++)
                for (size_t yoff = y * root, x = 0; x < root; x++)
                    deltas[x + yoff + zoff] = {delta_vals[x], delta_vals[y], delta_vals[z]};
    }
};

/*************************************************************************
  dim = 3
  interp = 1
  nbrRefineParticles = 27
*/
template<>
class Splitter<3, 1, 27> : public ASplitter<3, 1, 27>
{
public:
    Splitter()
    // : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    // constexpr static float delta_vals_[3]  = {{-.1, 0, .1}};
    // constexpr static float weight_vals_[3] = {{0.25, 0.125, 0.02725}};
};

/*************************************************************************
  dim = 3
  interp = 3
  nbrRefineParticles = 27
*/
template<>
class Splitter<3, 2, 27> : public ASplitter<3, 2, 27>
{
public:
    Splitter()
    // : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    // constexpr static float delta_vals_[3]  = {{-.1, 0, .1}};
    // constexpr static float weight_vals_[3] = {{0.25, 0.125, 0.02725}};
};


/*************************************************************************
  dim = 3
  interp = 3
  nbrRefineParticles = 27
*/
template<>
class Splitter<3, 3, 27> : public ASplitter<3, 3, 27>
{
public:
    Splitter()
    // : ASplitter{weight_vals_, delta_vals_}
    {
    }

protected:
    // constexpr static float delta_vals_[3]  = {{-.1, 0, .1}};
    // constexpr static float weight_vals_[3] = {{0.25, 0.125, 0.02725}};
};

} // namespace PHARE::amr
#endif /*PHARE_SPLIT_3D_H*/