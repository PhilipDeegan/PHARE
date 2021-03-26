
if(DEFINED CALIPER_ROOT)
  set(withCaliper ON)
endif()

if (withCaliper)
  if(DEFINED CALIPER_ROOT)
    find_package(caliper PATHS ${CALIPER_ROOT} REQUIRED)
    include_directories("${caliper_INCLUDE_DIR}")
  else()
    set(CALIPER_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/caliper)
    set(CALIPER_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/caliper)

    if (NOT EXISTS ${c})
      execute_process(
        COMMAND ${Git} clone https://github.com/LLNL/Caliper ${CALIPER_SRCDIR} -b master --recursive --depth 10
      )
    else()
      if(devMode)
        message("downloading latest Caliper updates")
        execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${CALIPER_SRCDIR})
      endif(devMode)
    endif()

    option(CALIPER_OPTION_PREFIX ON)
    set(CALIPER_OPTION_PREFIX ON)

    option(CALIPER_WITH_MPI ON)
    set(CALIPER_WITH_MPI ON)

    option(ENABLE_TESTS "Enable Samrai Test" OFF ) # disable CALIPER Test so that we can use the googletest pulled after

    add_subdirectory(${CALIPER_SRCDIR})
    unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE) # undoes what caliper does, so ctest can continue to work
    include_directories(${CMAKE_BINARY_DIR}/subprojects/caliper/include) # this is needed to find build-dir/include/SAMRAI/SAMRAI_config.h
    include_directories(${CALIPER_SRCDIR}/include)

  endif()

  add_definitions(-DPHARE_WITH_CALIPER=1)
  set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} caliper)
endif()
