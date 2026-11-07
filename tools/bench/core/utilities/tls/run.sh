#!/usr/bin/env bash
# thread_local access cost in a dlopened -fPIC module, see repo-root `next` (2026-11-03 item 8).
#
# Matrix: compiler {g++, clang++} x TLS dialect {gnu, gnu2} x visibility {default, hidden};
#  each matrix entry = one module with every access kind (tls_kind.cpp, one TU per kind).
#
#   ./run.sh syntax   # -fsyntax-only for every kind/flag combination, no objects
#   ./run.sh build    # objects + modules + driver into $OUT
#   ./run.sh sites    # per kind TLS relocation / call-site counts (needs build)
#   ./run.sh run      # time all modules (needs build), args after `run` go to the driver
#   ./run.sh all      # build + sites + run
#
# env: OUT (default build/tls_bench), CXXS (default "g++ clang++"), CXXFLAGS_EXTRA
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../../../.." && pwd)"
OUT="${OUT:-$ROOT/build/tls_bench}"
CXXS="${CXXS:-g++ clang++}"
DIALECTS="gnu gnu2"
VISIBILITIES="default hidden"
KINDS="0 1 2 3 4 5 6"
BASE_FLAGS="-std=c++20 -O3 -march=native -mtune=native -fPIC -Wall -Wextra -Werror ${CXXFLAGS_EXTRA:-}"

flags_for() { # dialect visibility
    echo "$BASE_FLAGS -mtls-dialect=$1 -fvisibility=$2"
}

mod_name() { # cxx dialect visibility
    echo "libtls_${1//+/x}_${2}_${3}"
}

do_syntax() {
    for cxx in $CXXS; do
        for d in $DIALECTS; do
            for v in $VISIBILITIES; do
                for k in $KINDS; do
                    # shellcheck disable=SC2046
                    $cxx $(flags_for "$d" "$v") -DTLS_KIND="$k" -fsyntax-only "$HERE/tls_kind.cpp"
                done
            done
        done
        $cxx $BASE_FLAGS -fsyntax-only "$HERE/tls_dummy.cpp"
        $cxx $BASE_FLAGS -fsyntax-only "$HERE/tls_bench.cpp"
        echo "syntax OK: $cxx"
    done
}

do_build() {
    mkdir -p "$OUT"
    for cxx in $CXXS; do
        for d in $DIALECTS; do
            for v in $VISIBILITIES; do
                local name objs=()
                name="$(mod_name "$cxx" "$d" "$v")"
                for k in $KINDS; do
                    local obj="$OUT/${name}_k$k.o"
                    # shellcheck disable=SC2046
                    $cxx $(flags_for "$d" "$v") -DTLS_KIND="$k" -c "$HERE/tls_kind.cpp" -o "$obj"
                    objs+=("$obj")
                done
                # shellcheck disable=SC2046
                $cxx $(flags_for "$d" "$v") -shared "${objs[@]}" -o "$OUT/$name.so"
            done
        done
    done
    local cxx0="${CXXS%% *}"
    $cxx0 $BASE_FLAGS -shared "$HERE/tls_dummy.cpp" -o "$OUT/libtls_dummy.so"
    $cxx0 $BASE_FLAGS -fno-pic -no-pie "$HERE/tls_bench.cpp" -o "$OUT/tls_bench" -ldl -pthread
    echo "built into $OUT"
}

do_sites() {
    # per object (= per kind): TLS relocation types (TLSGD/TLSLD = __tls_get_addr call site,
    #  GOTPC32_TLSDESC/TLSDESC_CALL = gnu2 descriptor call, GOTTPOFF = initial-exec),
    #  __tls_get_addr call count, TLS wrapper/init symbols, and the module's TLS segment size
    for so in "$OUT"/libtls_*_*_*.so; do
        local name
        name="$(basename "$so" .so)"
        echo "== $name  (PT_TLS memsz: $(readelf -lW "$so" | awk '$1=="TLS"{print $6}'))"
        for k in $KINDS; do
            local obj="$OUT/${name}_k$k.o"
            local relocs calls wrappers
            relocs="$(readelf -rW "$obj" | awk '$3 ~ /TLS|TPOFF|DTPOFF/ {print $3}' | sort | uniq -c \
                | awk '{printf "%s:%s ", $2, $1}')"
            calls="$(objdump -dr --no-show-raw-insn "$obj" | grep -c '__tls_get_addr' || true)"
            wrappers="$(nm -C "$obj" | grep -cE 'TLS (wrapper|init) function' || true)"
            printf "  k%s %-15s get_addr:%-3s wrappers:%-2s %s\n" "$k" \
                "$(grep -m1 -A1 "TLS_KIND == $k" "$HERE/tls_kind.cpp" | awk '/TLS_NAME/{print $3}')" \
                "$calls" "$wrappers" "$relocs"
        done
    done
}

do_run() {
    "$OUT/tls_bench" --dummy "$OUT/libtls_dummy.so" "$@" "$OUT"/libtls_*_*_*.so
}

cmd="${1:-all}"
shift || true
case "$cmd" in
    syntax) do_syntax ;;
    build) do_build ;;
    sites) do_sites ;;
    run) do_run "$@" ;;
    all) do_build && do_sites && do_run "$@" ;;
    *) echo "unknown command $cmd" >&2; exit 1 ;;
esac
