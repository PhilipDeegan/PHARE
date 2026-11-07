add_subdirectory(copy)
add_subdirectory(copy_overlap)
add_subdirectory(stream_pack)
add_subdirectory(refine)

function(_add_amr_particles_test src_name)

    add_executable(${src_name} ${src_name}.cpp)

    target_include_directories(${src_name} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      ${GTEST_INCLUDE_DIRS}
    )

    target_link_directories(${src_name} PRIVATE
      ${HDF5_LIBRARY_PATH}
    )

    target_link_libraries(${src_name} PRIVATE
      phare_amr
      ${GTEST_LIBS}
    )

    add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_amr_particles_test)

file(COPY "test_particles_schedules_inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

_add_amr_particles_test(test_particles_data)

# one binary per simulator permutation to prevent OOM when linking: several hierarchies
#  in one binary can exceed github actions runner memory (16GB), see tests/diagnostic
function(_add_particles_schedules_test dim layout L0 L1)
    set(src_name test_particles_schedules_${dim}d_${layout})

    add_executable(${src_name} test_particles_schedules.cpp)
    target_compile_definitions(${src_name} PRIVATE
      PHARE_TEST_DIM=${dim} PHARE_TEST_LAYOUT=${layout} PHARE_TEST_L0=${L0} PHARE_TEST_L1=${L1})

    target_include_directories(${src_name} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      ${GTEST_INCLUDE_DIRS}
    )

    target_link_directories(${src_name} PRIVATE
      ${HDF5_LIBRARY_PATH}
    )

    target_link_libraries(${src_name} PRIVATE
      phare_amr
      ${GTEST_LIBS}
    )

    add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_particles_schedules_test)

#                             dim layout    L0 L1
_add_particles_schedules_test(1   AoSMapped 1  1)
_add_particles_schedules_test(2   AoSMapped 0  1)
if(withMkn OR withGPU)
  _add_particles_schedules_test(1 AoSPCTS   1  0)
endif()

