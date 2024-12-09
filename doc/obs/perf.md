## Performance analyses

### via perf

See: https://man7.org/linux/man-pages/man1/perf-record.1.html

```shell
# build with the following
CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer"
CMAKE_BUILD_TYPE="RelWithDebInfo"

# build...

perf record -Tga -F 1000 ./build/src/phare/phare-exe path/to/python/script.py
ls -lah perf.data # check file is new/exists
perf script report flamegraph
```

### via perf in vtune

installation see
https://www.intel.com/content/www/us/en/docs/vtune-profiler/installation-guide/2023-0/package-managers.html


```shell
# see above, but rename perf.data to data.perf (vtune only reads .perf files)
```

### via scalasca / scorep

```shell
# install
sudo dnf install scorep-openmpi scalasca-openmpi libunwind-devel binutils-devel elfutils-devel
```

```shell
# usage
CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer"
CMAKE_BUILD_TYPE="RelWithDebInfo"
export CXX=scorep-mpicxx CC=scorep-mpicc FC=scorep-gfortran

# build...

# run
# see https://vampir.eu/public/files/pdf/spcheatsheet_a4.pdf
export SCOREP_EXPERIMENT_DIRECTORY=scorep_run_trace
export SCOREP_TOTAL_MEMORY=500M # default is 16MB which will fail
scalasca -analyze mpirun -n 8 ./build/src/phare/phare-exe path/to/python/script.py

# view
scalasca -examine scorep_run_trace
```
