// compiled with the GPU flags the configurator would give PHARE (see config.py): proves the
//  compiler can target the device, the kernel runs and thrust (needed by PHARE's device
//  sorts) works, exit code 0 on success

#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#if PHARE_PROBE_HIP
#include <hip/hip_runtime.h>
#define PROBE_GPU(x) hip##x
#else
#include <cuda_runtime.h>
#define PROBE_GPU(x) cuda##x
#endif

__global__ void phare_gpu_check(int* const v)
{
    int a[3] = {42, 1, 2};
    thrust::sort(thrust::seq, a, a + 3);
    *v = a[2];
}

int main()
{
    int* v = nullptr;
    if (PROBE_GPU(MallocManaged)(&v, sizeof(int)) != PROBE_GPU(Success))
        return 1;
    *v = 0;
    phare_gpu_check<<<1, 1>>>(v);
    if (PROBE_GPU(DeviceSynchronize)() != PROBE_GPU(Success))
        return 2;
    int const ok = *v == 42 ? 0 : 3;
    PROBE_GPU(Free)(v);
    return ok;
}
