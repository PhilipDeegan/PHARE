cmake_minimum_required (VERSION 3.20.1)

project(test-updater)

add_executable(${PROJECT_NAME} test_updater.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )
target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  phare_simulator
  ${GTEST_LIBS})
add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_updater_fixtures test_fixtures.cpp)
target_include_directories(test_updater_fixtures PRIVATE
  ${GTEST_INCLUDE_DIRS}
  )
target_link_libraries(test_updater_fixtures PRIVATE
  phare_core
  phare_simulator
  ${GTEST_LIBS})
add_no_mpi_phare_test(test_updater_fixtures ${CMAKE_CURRENT_BINARY_DIR})

