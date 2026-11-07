## PHARE Configurator

This tool will investigate your system and setup PHARE to compile in the most optimal way for your hardware.

Such that if GPUs are available (and working), PHARE will be compiled to work with those provided there is a compatible compiler available.

Run at CMake configure time with `-Dphare_configurator=ON`, output in `${CMAKE_BINARY_DIR}/phare_configurator`:
- MPI / HDF5 versions (main.cpp), the generated headers in `${CMAKE_BINARY_DIR}/src/core/def`
- GPU (config.py): CUDA then HIP toolkit, device arch from a host only runtime query (gpu_probe.cpp),
  then a test kernel built with the flags PHARE would use and run (gpu_check.cpp)
- `local.cmake`, included by `res/cmake/def.cmake`: `PHARE_GPU_VENDOR`, `PHARE_GPU_ARCH`,
  `PHARE_GPU_CXX_FLAGS`, `PHARE_GPU_LINK_FLAGS` (empty if no usable GPU + compiler)

The GPU compiler is clang (`-x cuda` / `-x hip`), nvcc is not supported: with any other
`CMAKE_CXX_COMPILER` no GPU is configured. The flags only apply with `-DwithGPU=ON` (opt-in,
res/cmake/gpu.cmake), which fails the configure if no usable GPU + compiler was found.

Overrides: `-DPHARE_GPU=auto|off|cuda|hip`, `-DPHARE_GPU_ARCH=sm_90|gfx90a|...` (skips the device
query, e.g. on a build node without a GPU: the test kernel is then compiled but not run).
