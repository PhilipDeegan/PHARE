#!/usr/bin/env bash

set -ex
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD

#time (
  date

  # FILE="-M tests/core/data/particles/test_particles.cpp"
  # FILE="-M tests/amr/data/particles/stream_pack/test_main.cpp"
  # FILE="-M tests/core/data/particles/sorting/test_particle_sorting.cpp"
  # FILE="-M tests/core/data/particles/test_edge_bisection_mapper.cpp"
  # FILE="-M tests/core/data/particles/test_bisection_range_mapper.cpp"
  FILE="-M tests/core/utilities/range/test_range.cpp"
  FILE="-M tests/core/numerics/interpolator/test_main.cpp"
  FILE="-M tests/core/numerics/pusher/test_pusher.cpp"
  ARGS="${FILE} -l -pthread " #-x res/mkn/clang_asan"

  # mkn clean build run -p test_amr $ARGS -WgO 9 || true
  mkn clean build -p test_core dbg $ARGS -WgO  $@ || true
  date
# ) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
#) 1> $CWD/.mkn.sh.out 2> $CWD/.mkn.sh.err >&2



