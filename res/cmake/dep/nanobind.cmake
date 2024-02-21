
SET(PYBIND_MIN_VERSION "2.5.0")

function(get_nanobind)

  message("downloading subproject nanobind")
  set(NANOBIND_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/nanobind)

  if (NOT EXISTS ${NANOBIND_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/wjakob/nanobind ${NANOBIND_SRCDIR} --depth 1 -b master
    )
  else()
    if(devMode)
      message("downloading latest nanobind updates")
      execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${NANOBIND_SRCDIR})
    endif(devMode)
  endif()

  add_subdirectory(${NANOBIND_SRCDIR})

endfunction(get_nanobind)

if (forceGetNanobind)
  get_nanobind()
else()

  find_package(nanobind ${PYBIND_MIN_VERSION} CONFIG QUIET)

  if (NOT nanobind_FOUND)
    get_nanobind()
  endif()

endif(forceGetNanobind)
