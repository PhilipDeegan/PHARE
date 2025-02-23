
#include "evalOnBox.hpp"


#include "mkn/avx/span.hpp"
#include "mkn/avx/vector.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    using Array_t = mkn::avx::Vector_t<T>;

    GridLayout_t layout{C - 1};
    auto const local_box
        = PHARE::core::box_from_zero_to_upper(*layout.AMRBox().shape().as_unsigned());
    std::array<std::uint32_t, dim> const dims = local_box.shape();
    auto const size                           = PHARE::core::product(dims);

    Array_t const rho(size, 1);
    Array_t tmp(size, 0);

    auto const v0 = mkn::avx::make_span(rho);
    auto v1       = mkn::avx::make_span(tmp);

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            v1 += v0;
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
