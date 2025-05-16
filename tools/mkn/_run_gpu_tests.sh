# expects defined functions run and run_mpi # i.e. see cuda.sh


run tests/core/data/test_vector.cpp
run tests/core/data/particles/sorting/test_gpu_sorting_mkn.cpp
run tests/core/numerics/pusher/test_gpu_pusher.cpp
run tests/core/numerics/faraday/test_gpu_main.cpp
run tests/core/numerics/ampere/test_gpu_main.cpp
run tests/core/numerics/ohm/test_gpu_main.cpp

run_mpi tests/core/data/test_gpu_mpi_main.cpp 3
