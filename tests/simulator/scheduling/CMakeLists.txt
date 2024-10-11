cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-scheduling)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

function(_add_scheduling_test src_name)

  add_executable(${src_name} ${src_name}.cpp)

  target_include_directories(${src_name} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_directories(${src_name} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )

  target_link_libraries(${src_name} PRIVATE
    phare_simulator pybind11::embed phare_initializer
    ${HDF5_LIBRARIES})

  add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_scheduling_test)

if(HighFive)
  if(testMPI)
    _add_scheduling_test(test_simulator_scheduling )
  endif(testMPI)
endif()
