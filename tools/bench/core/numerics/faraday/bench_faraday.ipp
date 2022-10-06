/*
#include "bench/core/numerics/faraday/bench_faraday.ipp"
*/
#include "bench/core/bench.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/faraday/faraday_avx.hpp"

#include "mkn/kul/time.hpp"

constexpr std::size_t N       = 1;
constexpr std::size_t dim     = 3;
constexpr std::size_t interp  = 1;
constexpr std::uint32_t cells = 300;
constexpr double dt           = .001;

using PHARE_Types  = PHARE::core::PHARE_Types<dim, interp>;
using GridLayout_t = typename PHARE_Types::GridLayout_t;
using PHARE::core::Component;
using PHARE::core::Direction;