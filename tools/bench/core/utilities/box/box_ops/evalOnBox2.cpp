
#include "evalOnBox.hpp"
#include "core/utilities/box/box_spans.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    GridLayout_t layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    Grid_t tmp{"tmp", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            PHARE::core::for_each_point(layout.domainBoxFor(rho),
                                        [&](auto const& bix) { tmp(bix) += rho(*bix); });
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
