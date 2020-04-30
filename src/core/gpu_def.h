#ifndef PHARE_CORE_GPU_DEF_H
#define PHARE_CORE_GPU_DEF_H

#if defined(PHARE_WITH_HIP)

#define PHARE_GPU_DEV_ONLY __device__
#define PHARE_GPU_HST_ONLY __device__
#define PHARE_GPU_FUNC PHARE_GPU_HST_ONLY PHARE_GPU_DEV_ONLY

#else

#define PHARE_GPU_DEV_ONLY
#define PHARE_GPU_HST_ONLY
#define PHARE_GPU_FUNC

#endif

#endif /*PHARE_CORE_GPU_DEF_H*/
