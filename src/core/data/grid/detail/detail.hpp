#ifndef PHARE_CORE_GRID_DETAIL_DETAIL_HPP
#define PHARE_CORE_GRID_DETAIL_DETAIL_HPP

#if PHARE_HAVE_RAJA
#include "core/data/grid/detail/raja.hpp"
#endif

#if PHARE_HAVE_MKN_GPU
#include "core/data/grid/detail/mkn_gpu.hpp"
#endif


#endif /* PHARE_CORE_GRID_DETAIL_DETAIL_HPP */
