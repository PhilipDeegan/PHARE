
if (withPhlop)

  # pinned by commit: the python package (pyproject.toml "timing") must read what these headers write
  # temporary: PhilipDeegan/phlop next (#64-#76 + scope_compare runs x ranks)
  if(NOT DEFINED PHARE_PHLOP_VERSION)
    SET(PHARE_PHLOP_VERSION "3e0a580e299565b77face24ea66637d63e97da61")
  endif()

  set(PHLOP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/phlop)
  if (NOT EXISTS ${PHLOP_SRCDIR})
    execute_process(COMMAND ${Git} init -q ${PHLOP_SRCDIR})
    execute_process(COMMAND ${Git} remote add origin https://github.com/PhilipDeegan/phlop
                    WORKING_DIRECTORY ${PHLOP_SRCDIR})
  endif()

  execute_process(COMMAND ${Git} rev-parse HEAD WORKING_DIRECTORY ${PHLOP_SRCDIR}
                  OUTPUT_VARIABLE PHLOP_HEAD OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if (NOT PHLOP_HEAD STREQUAL PHARE_PHLOP_VERSION)
    message("fetching phlop ${PHARE_PHLOP_VERSION}")
    execute_process(COMMAND ${Git} fetch -q --depth 1 origin ${PHARE_PHLOP_VERSION}
                    WORKING_DIRECTORY ${PHLOP_SRCDIR} COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND ${Git} checkout -q FETCH_HEAD
                    WORKING_DIRECTORY ${PHLOP_SRCDIR} COMMAND_ERROR_IS_FATAL ANY)
  endif()

  include_directories(${PHLOP_SRCDIR}/inc)

endif(withPhlop)
