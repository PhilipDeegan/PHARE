#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e
. $CWD/tools/mkn_func.sh
CARGS=${CARGS:-""}

mkdir -p .phare/timings/tests/amr/data/field

# FILE="tests/core/numerics/ion_updater/test_updater_pp_main.cpp"
# FILE="tests/core/numerics/faraday/test_tileset_faraday.cpp"
# FILE="tests/diagnostic/test_diagnostics_1d_cmp.cpp"
# FILE="tests/core/numerics/ampere/test_tileset_ampere.cpp"
# FILE="tests/diagnostic/test-diagnostics_2d.cpp"
# FILE="tests/core/numerics/test_field_evolution.cpp"
# FILE="tests/amr/data/particles/test_particles_schedules.cpp"
# FILE="tests/amr/data/particles/refine/test_particles_data_split.cpp"
# FILE="tests/diagnostic/test-diagnostics_1d.cpp"
FILE="tests/core/data/particles/test_particles_partitioning.cpp"
FILE="tests/core/data/particles/test_particles_AoSTS.cpp"
FILE="tests/core/data/grid/test_grid_tile_set.cpp"
FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
FILE="tools/bench/amr/data/field/schedules/bench_field_schedules.cpp"
FILE="tests/core/data/utilities/box/test_box_span.cpp"
FILE="tests/core/data/field/test_field_overlaps.cpp"
FILE="tests/core/data/particles/test_particles_construction.cpp"
FILE="tests/amr/data/field/test_L0_fields_schedules.cpp"
FILE="tests/amr/data/particles/copy/test_particledata_copyNd.cpp"
FILE="tests/amr/data/field/test_fields_schedules.cpp"
FILE="tests/amr/data/field/test_L0_fields_schedules.cpp"
FILE="tests/amr/data/particles/test_particles_data.cpp"
FILE="tests/core/numerics/ion_updater/test_updater.cpp"
FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
TEST=" -M ${FILE} "
CARGS="${CARGS} $(mkn_get_opts)"

set -x
(
  [[ $FILE == tests/core/* ]] && (
    [ ! -f "bin/core/libphare_core.a" ] && (
        mkn clean build -KOgqtp test_core -d google.test,+ -a "-fPIC" ${CARGS}
        mkn clean build -KOgqtp core -a "-fPIC" ${CARGS}
    )
    mkn -p test_core ${TEST} ${CARGS} "$@"
  ) || true

  [[ $FILE != tests/core/* ]] && (
    [ ! -f "bin/amr/libphare_amr.a" ] && (
        mkn clean build -KOgqtdp test_amr -a "-fPIC" ${CARGS}
    )
    mkn -qp test_diagnostics ${TEST} ${CARGS} "$@"
    echo "ok"
  ) || true

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

#
TIME_FILE=".phare/timings/tests/amr/data/field/test_L0_fields_schedules.cpp.txt"
[ -f "${TIME_FILE}" ] && \
  python3 tools/python3/phloping.py print_scope_timings -f "${TIME_FILE}" | tee times.txt
#


# stop if not CI
[ -z "$CI" ] && exit 0 # continue for more

(
  set -xe
  mkn -tp core_tests ${CARGS} clean build test run -gOw mkn.gpu
  mkn -tp more_tests ${CARGS} clean build test run -gOw mkn.gpu
)
