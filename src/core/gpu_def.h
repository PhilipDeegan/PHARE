#ifndef PHARE_CORE_GPU_DEF_H
#define PHARE_CORE_GPU_DEF_H

#if defined(PHARE_WITH_HIP)

#define PHARE_GPU_DEVICE __device__
#define PHARE_GPU_HST __host__
#define _PHARE_FN_DECORATORS_ PHARE_GPU_HST PHARE_GPU_DEVICE

#else

#define PHARE_GPU_DEVICE
#define PHARE_GPU_HST
#define _PHARE_FN_DECORATORS_

#endif

#endif /*PHARE_CORE_GPU_DEF_H*/
