set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD # move to project root

MAIN="-M tools/bench/core/numerics/ion_updater/bench_ion_updater.cpp"
time (
  KLOG=0 mkn clean build ${MAIN} -p bench run -Og 0
) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )
