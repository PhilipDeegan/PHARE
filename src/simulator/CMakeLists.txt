cmake_minimum_required (VERSION 3.20.1)
project(phare_simulator)

set( SOURCES_INC
    simulator.h
    phare_types.h
   )


add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/phare/simulator>)

target_link_libraries(${PROJECT_NAME} INTERFACE
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
)
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${PHARE_LINK_FLAGS}")

