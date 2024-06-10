
#  TO activate phlop with PHARE
#   1. configure cmake -DwithPhlop=ON
#   2. export PHARE_SCOPE_TIMING=1
#     for function scope timings
##

if (withPhlop)

  set(PHLOP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/phlop)

  if (NOT EXISTS ${PHLOP_SRCDIR})
    execute_process(
      COMMAND ${Git} clone https://github.com/PhilipDeegan/phlop ${PHLOP_SRCDIR} -b master --recursive --depth 1
    )
  else()
    if(devMode)
      message("downloading latest phlop updates")
      execute_process(COMMAND ${Git} pull origin master WORKING_DIRECTORY ${PHLOP_SRCDIR})
    endif(devMode)
  endif()

  include_directories(${PHLOP_SRCDIR}/inc)

endif(withPhlop)

