#  Darshan (ANL) I/O characterization - https://github.com/darshan-hpc/darshan
#
#  TO activate darshan with PHARE
#   1. configure cmake -DwithDarshan=ON (optionally -DHDF5_ROOT=/path for the HDF5 module)
#   2. build, darshan is installed to ${CMAKE_BINARY_DIR}/subprojects/darshan/install
#   3. run with LD_PRELOAD, e.g.
#       mpirun -x LD_PRELOAD=${PHARE_DARSHAN_LIB} -n 4 python3 job.py
#   4. logs go to $DARSHAN_LOG_DIR if set, else the current working directory
#       ${DARSHAN_INSTALL}/bin/darshan-parser <logfile>
##

if (withDarshan)

  if(NOT DEFINED PHARE_DARSHAN_VERSION)
    set(PHARE_DARSHAN_VERSION "darshan-3.4.7")
  endif()

  set(DARSHAN_SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/darshan)
  set(DARSHAN_BIN ${CMAKE_CURRENT_BINARY_DIR}/subprojects/darshan)
  set(DARSHAN_INSTALL ${DARSHAN_BIN}/install)

  phare_github_get_or_update(darshan ${DARSHAN_SRCDIR} darshan-hpc/darshan ${PHARE_DARSHAN_VERSION})

  set(DARSHAN_CONFIGURE_ARGS
    --prefix=${DARSHAN_INSTALL}
    --with-log-path-by-env=DARSHAN_LOG_DIR,PWD
    --with-jobid-env=NONE
    CC=${MPI_C_COMPILER}
  )

  # darshan expects an install prefix with bin/h5dump and lib/, so only with HDF5_ROOT
  if(DEFINED HDF5_ROOT AND HDF5_IS_PARALLEL)
    list(APPEND DARSHAN_CONFIGURE_ARGS --with-hdf5=${HDF5_ROOT})
    message("darshan HDF5 module enabled with HDF5_ROOT: ${HDF5_ROOT}")
  else()
    message("darshan HDF5 module disabled, requires -DHDF5_ROOT for parallel HDF5")
  endif()

  include(ExternalProject)
  ExternalProject_Add(darshan
    SOURCE_DIR ${DARSHAN_SRCDIR}
    BINARY_DIR ${DARSHAN_BIN}/build
    INSTALL_DIR ${DARSHAN_INSTALL}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir ${DARSHAN_SRCDIR} ./prepare.sh
              COMMAND ${DARSHAN_SRCDIR}/configure ${DARSHAN_CONFIGURE_ARGS}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )

  set(PHARE_DARSHAN_LIB ${DARSHAN_INSTALL}/lib/libdarshan.so)
  message("darshan LD_PRELOAD library: ${PHARE_DARSHAN_LIB}")

endif(withDarshan)
