

if (withPhlop)

  set(PHLOP_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/phlop)
  if(NOT DEFINED PHARE_PHLOP_BRANCH)
    SET(PHARE_PHLOP_BRANCH "master")
  endif()

  # devMode pulls ${PHARE_PHLOP_BRANCH} into an existing checkout, only do so when that
  #  branch is checked out, so a local working branch (eg next) is never merged into
  set(PHLOP_CURRENT_BRANCH "")
  if (EXISTS ${PHLOP_SRCDIR})
    execute_process(COMMAND ${Git} branch --show-current WORKING_DIRECTORY ${PHLOP_SRCDIR}
                    OUTPUT_VARIABLE PHLOP_CURRENT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  endif()

  if (EXISTS ${PHLOP_SRCDIR} AND NOT PHLOP_CURRENT_BRANCH STREQUAL PHARE_PHLOP_BRANCH)
    message("phlop: on '${PHLOP_CURRENT_BRANCH}', not '${PHARE_PHLOP_BRANCH}', not updating")
  else()
    phare_github_get_or_update(phlop ${PHLOP_SRCDIR} PhilipDeegan/phlop ${PHARE_PHLOP_BRANCH})
  endif()

  include_directories(${PHLOP_SRCDIR}/inc)

endif(withPhlop)
