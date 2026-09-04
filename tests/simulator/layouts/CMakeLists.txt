cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-layouts)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive AND withMkn)
  ## These tests use dump diagnostics so require HighFive, and compare
  ## against mkn-only particle layouts so require withMkn!

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 1)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} layouts-1d test_layouts_1d.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 2)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} layouts-2d test_layouts_2d.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 3)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} layouts-3d test_layouts_3d.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()

endif()
