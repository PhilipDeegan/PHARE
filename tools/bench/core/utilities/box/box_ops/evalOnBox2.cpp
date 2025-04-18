
#include "evalOnBox.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    GridLayout_t layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    Grid_t tmp{"tmp", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            for (auto const& bix : layout.domainBoxFor(rho))
                tmp(bix) += rho(*bix);
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
