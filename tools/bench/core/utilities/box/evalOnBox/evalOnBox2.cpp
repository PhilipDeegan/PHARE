
#include "evalOnBox.hpp"
#include "core/utilities/box/box_spans.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    GridLayout_t layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};

    T x = 0;
    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            PHARE::core::for_each_point(layout.domainBoxFor(rho),
                                        [&](auto const& bix) { x += rho(*bix); });
    }

    if (x != C * C * C * N)
        throw std::runtime_error("no");
}
