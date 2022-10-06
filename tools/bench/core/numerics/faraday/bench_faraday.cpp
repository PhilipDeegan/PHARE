#include "bench/core/numerics/faraday/bench_faraday.ipp"

int main(int argc, char** argv)
{
    KOUT(NON) << __FILE__;

    auto layout = PHARE::core::bench::getLayout<dim, interp>(cells);

    PHARE::core::bench::Electromag EM{layout}, EMNew{layout};
    PHARE::core::Faraday<GridLayout_t> faraday;
    auto __ = PHARE::core::SetLayout(&layout, faraday);

    auto const s = mkn::kul::Now::MILLIS();
    for (std::size_t i = 0; i < N; ++i)
        faraday(EM.B, EM.E, EMNew.B, dt);
    KOUT(NON) << "RUN TIME: " << (mkn::kul::Now::MILLIS() - s) << " ms";
}
