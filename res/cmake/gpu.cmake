# -DwithGPU=ON (opt-in): compile and link every PHARE target for the GPU found by the
#  configurator (tools/config/config.py, PHARE_GPU_* from its local.cmake). Clang only (-x cuda /
#  -x hip), no nvcc. Included after dep.cmake, so third party subprojects (SAMRAI, highfive,
#  gtest, ...) are not compiled as CUDA/HIP: add_compile_options only reaches directories added
#  after this

if(withGPU)
  if(NOT phare_configurator)
    message(FATAL_ERROR "withGPU needs -Dphare_configurator=ON (GPU detection)")
  endif()
  if(NOT PHARE_GPU_VENDOR)
    message(FATAL_ERROR "withGPU: no usable GPU + compiler, configurator says: ${PHARE_GPU_STATUS}"
                        " (see -DPHARE_GPU / -DPHARE_GPU_ARCH, tools/config/README.md)")
  endif()

  if("${CMAKE_CXX_FLAGS}" MATCHES "MKN_GPU_CPU")
    message(FATAL_ERROR "withGPU: CMAKE_CXX_FLAGS has MKN_GPU_CPU, which forces mkn.gpu's CPU backend")
  endif()
  if("${CMAKE_CXX_FLAGS}" MATCHES "PHARE_HAVE_THRUST=0")
    message(FATAL_ERROR "withGPU: CMAKE_CXX_FLAGS has PHARE_HAVE_THRUST=0, device sorts need thrust")
  endif()

  message(STATUS "PHARE GPU build: ${PHARE_GPU_STATUS}")
  add_compile_options(${PHARE_GPU_CXX_FLAGS})
  add_link_options(${PHARE_GPU_LINK_FLAGS})
  # clang warns for nvcc (_PHARE_ALL_FN_ after a lambda's parameter list), see devMode Werror
  check_cxx_compiler_flag(-Wno-cuda-compat PHARE_HAS_WNO_CUDA_COMPAT)
  if(PHARE_HAS_WNO_CUDA_COMPAT)
    add_compile_options(-Wno-cuda-compat)
  endif()
elseif(withMkn)
  # mkn.gpu picks its backend from the headers it finds (e.g. CUDA if installed): without
  #  withGPU nothing is compiled for the device, so force its CPU backend (thrust is off
  #  already: core/def/detail/thrust.hpp needs PHARE_HAVE_GPU)
  add_compile_definitions(MKN_GPU_CPU=1)
endif(withGPU)
