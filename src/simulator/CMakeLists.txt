cmake_minimum_required(VERSION 3.9)
project(phare_simulator)

set( SOURCES_INC
    simulator.h
    phare_types.h
   )

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE
  phare_initializer
  phare_amr # for mpicc
  phare_diagnostic
  )

