cmake_minimum_required (VERSION 3.20.1)

project(test-pharesee)


add_python3_test(test-pharesee-geometry_1d test_geometry.py ${PROJECT_SOURCE_DIR})
add_python3_test(test-pharesee-geometry_2d test_geometry_2d.py ${PROJECT_SOURCE_DIR})
add_python3_test(test-pharesee-hierarchy test_hierarchy.py ${PROJECT_SOURCE_DIR})
#phare_mpi_python3_exec(11 2 test-pharesee-hierarchy test_hierarchy.py ${PROJECT_SOURCE_DIR})


