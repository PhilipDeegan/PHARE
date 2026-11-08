// host only: queries the GPU runtime, compiled by the plain C++ compiler (see config.py)
//  one line per device, tab separated key=value

#include <cstdio>

#if PHARE_PROBE_HIP
#include <hip/hip_runtime_api.h>
#define PROBE_GPU(x) hip##x
using DeviceProp_t = hipDeviceProp_t;
#else
#include <cuda_runtime_api.h>
#define PROBE_GPU(x) cuda##x
using DeviceProp_t = cudaDeviceProp;
#endif

int main()
{
    int n = 0;
    if (PROBE_GPU(GetDeviceCount)(&n) != PROBE_GPU(Success) or n == 0)
        return 1;

    for (int i = 0; i < n; ++i)
    {
        DeviceProp_t p;
        if (PROBE_GPU(GetDeviceProperties)(&p, i) != PROBE_GPU(Success))
            return 2;

#if PHARE_PROBE_HIP
        std::printf("device=%d\tname=%s\tarch=%s\twarp=%d\tmax_threads=%d\tmanaged=%d\n", i, p.name,
                    p.gcnArchName, p.warpSize, p.maxThreadsPerBlock, p.managedMemory);
#else
        std::printf("device=%d\tname=%s\tarch=sm_%d%d\twarp=%d\tmax_threads=%d\tsmem_optin=%zu\t"
                    "managed=%d\n",
                    i, p.name, p.major, p.minor, p.warpSize, p.maxThreadsPerBlock,
                    p.sharedMemPerBlockOptin, p.managedMemory);
#endif
    }
    return 0;
}
