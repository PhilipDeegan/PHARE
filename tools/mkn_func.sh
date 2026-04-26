

mkn_rocm_available(){
  V=0
  # [ -d /opt/rocm/bin ] && V=$((V + 1))
  # (( V == 0 )) && which hipcc 2>&1 > /dev/null && V=1
  echo $V
}


mkn_cuda_available(){
  V=0
  # which clang 2>&1 > /dev/null && clang -v 2>&1 | grep -q "Found CUDA" && V=1
  echo $V
}


mkn_get_opts () {
  set -e
  CUDA="$(mkn_cuda_available)"
  ROCM=$(mkn_rocm_available)
  (( $ROCM > 0 )) && (( $CUDA > 0 )) && echo "ERROR: BOTH ROCM AND CUDA FOUND!" && exit 1
  ARGS="" # enable gpu if support is considered available
  (( $ROCM > 0 )) && ARGS="res/mkn/hip_mpi"
  (( $CUDA > 0 )) && ARGS+="res/mkn/clang_cuda"
  [ -n "$ARGS" ] && ARGS="-P mkn.base=gpu_ -x $ARGS"
  [ -z "$ARGS" ] && ARGS="-x res/mkn/mpi" # default
  echo "-w mkn.gpu $ARGS"
}
