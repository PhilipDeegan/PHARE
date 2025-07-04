


mkn_get_opts () {
  set -e
  ARGS="" # enable gpu if support is considered available
  # [ -d /opt/rocm/bin ] && ARGS="res/mkn/hip_mpi"
  # [ -z "$ARGS" ] && which clang 2>&1 > /dev/null && clang -v 2>&1 | grep -q "Found CUDA" && ARGS="res/mkn/clang_cuda"
  [ -n "$ARGS" ] && ARGS="-P mkn.base=gpu_ -x $ARGS"
  [ -z "$ARGS" ] && ARGS="-x res/mkn/mpi" # default
  echo "$ARGS"
}