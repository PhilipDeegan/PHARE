

if(NOT DEFINED PHARE_CPPDICT_VERSION)
  SET(PHARE_CPPDICT_VERSION "superserial")
endif()

set(CPPDICT_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/cppdict)
phare_github_get_or_update(cppdict ${CPPDICT_SRCDIR} PhilipDeegan/cppdict ${PHARE_CPPDICT_VERSION})
