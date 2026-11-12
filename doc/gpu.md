# GPU (tiling_gpu): status, findings, todo

Build: RAJA + Umpire + SAMRAI stack (clang CUDA, sm_120), `./tools/cmake.sh` (git-ignored, `MODULES`
list there). The stack is `/mkn/r/llnl/{umpire,raja,samrai}/phare` (each `mkn.sh`: `ENABLE_CUDA`,
`BLT_ENABLE_CLANG_CUDA`); `tools/cmake.sh` points at it with `RAJA_UMPIRE=1`. Layout tests:
`PYTHONPATH=$PWD:$PWD/build:$PWD/pyphare python3 -m tests.simulator.layouts.test_layouts_{1,2,3}d`.
Core tests and benchmarks iterate faster with mkn (`./mkn_cuda.sh`, `mkn ... -p bench`), but mkn's
`samrai.ver: master` is `/mkn/r/llnl/samrai/master` (no RAJA/CUDA), see Open > Build.

Building in `build/` directly: `ninja -C build -j4 <targets>`. pyphare builds a missing module on
import (`pyphare/pyphare/cpp/cmake.py`, `-j ${THREADS:-4}`): after a reconfigure, rebuild the modules
a test imports before running it.

## Where things stand (2026-11-12)

Layout tests (incl. SoAVXChunked CPU + GPU_UNIFIED): 2d 25/25, 3d 24/25 (AoSTS GPU L2 level ghost
bug, see Open). Most of 2026-11-11/12 below is uncommitted.

3D harris (AoSTS `GPU_UNIFIED`, 40x80x20 cells, 100 ppc, 20 steps): 2.25 -> 0.88 s/step (2026-11-10:
block per tile field kernels + device realloc 1.27, device `average_` + `sync_inner_ghosts` 0.88).

2d harris, 100x100 cells, 100 ppc, 3 levels, 1 coarse step, 1 rank, 1 thread (seconds):

| scope | AoSMapped CPU | SoAVXChunked CPU | SoAVXChunked GPU 11-12 am | SoAVXChunked GPU now |
|---|---|---|---|---|
| `Simulator::advance` | 2.57 | 4.72 | 21.95 | ~7.3-8.1 |
| `SolverPPC::moveIons` | 2.12 | 3.83 | 13.07 (host deposit) | 0.35 |
| `SolverPPC::average_` | 0.005 | 0.003 | 0.81 (host) | 0.023 |
| field `*::schedules` (pred/corr/average) | ~0.05 | ~0.05 | ~2.5 | ~2.5 |
| `regrid` / `firstStep` | 0.17 / 0.06 | 0.16 / 0.07 | 1.37 / 0.99 | 1.38 / 0.98 |

What's left on GPU is host side work on managed data, mostly in the messenger: see Open > UM faults.

## Done

### Field solvers and pointwise ops: one launch per level, tiled or not
- `Ampere`, `Faraday`, `Ohm` are point operations (`core/numerics/{ampere,faraday,ohm}/*.hpp`),
  `core::*Transformer` (`*_transformer.hpp`) evaluate them: CPU untiled per patch
  (`layout.evalOn*Box`), CPU tiled per tile (`on_tile`, pool task)
- `core/utilities/kernels/region_chunks.hpp`: `run_region_chunks(items, fn)`, items carry
  `boxes[comp]`, each (item, comp) box split evenly into ceil(size / target) chunks of
  ceil(size / n) points (no mostly idle remainder chunk), a block per chunk, one launch. Target
  from the device: ~8 blocks per SM within [1024, 16384] points (`gpu::sm_count()`)
- `untiled_view()` (`field_tiles.hpp`): untiled (vec)fields as `basic::` views, as `tile_at` gives
  tiles, so an untiled patch is one region and a tile set one region per tile
- Ampere: `AmpereTransformer::gpu_items` / `run_gpu`, `AmpereLevelTransformer::gpu_level_`
  collects every patch's regions, one launch, then (tiled) each patch's tile ghost sync. Faraday /
  Ohm GPU still: tiled one `FieldTileLauncher` kernel per patch, untiled `evalOnBox` per component
  per patch (synchronous launches). Their kernels are cheap (pred/corr kernels 0.01-0.1 s total),
  the time is in the schedules
- `average`: `core/utilities/algorithm.hpp` `average_jobs` (buffers: untiled field one, tile set
  one per tile, 4096 point chunks), `run_average_jobs<alloc_mode>` (device launch or host loop),
  `average_job_list` (reused per thread). `core::average` uses them for GPU_UNIFIED;
  `SolverPPC::average_` collects E and B of every patch, one launch per level.
  `average_on_device` (field_tiles_gpu.hpp) removed
- `sync_inner_ghosts_on_device` (tiled only), `FieldTileSet::for_each_inner_ghost_overlap`
- tests: `tests/core/numerics/{ampere,faraday,ohm}/test_tileset_*.cpp` (`device_ops` uses
  `core::average`)

### Ion update on the device for untiled GPU_UNIFIED (SoAVXChunked)
- `IonUpdaterFor_t`: GPU_UNIFIED, tiled or not -> `mkn_xyz::ParallelIonUpdater` (deposit with
  `particleToMesh` on each patch's stream, moments zeroed on the streams). Untiled GPU used to go
  through `FlatParallelIonUpdater`, which deposited on the host (5.1 s of 13 s moveIons)
- `ParallelIonUpdater::post_move_sync_` re-sorts SoAVXChunked (`on_moved`)
- level 2 scopes in `FlatParallelIonUpdater` phases and `SoAVXChunked::sort`

### SAMRAI kernel fusion for field copies
- SAMRAI schedules give fuseable transactions (`RefineCopyTransaction`, `CoarsenCopyTransaction`)
  `StagedKernelFusers`, and call `copyFuseable` / `packStreamFuseable` / `unpackStreamFuseable`,
  but never launch the fusers: that is a `tbox::ScheduleOpsStrategy`'s job
- `amr/messengers/fuser_ops_strategy.hpp`: `amr::SetEqualFuseable` (its `span(dst, src, n)` is
  enqueued in stage 0 with `StagedKernelFusers::enqueue`, only between a strategy pre/post hook on
  this thread, else copied at once), `amr::FuserOpsStrategy` (pre hooks open the window, post hooks
  launch if `enqueued()`). Must use the staged `enqueue`: `getFuser(0)->enqueue` doesn't flag
  `enqueued()`, so nothing launched and the copies ran at someone else's launch (illegal address)
- core hook only: an operator with a static `span(dst, src, n)` gets whole spans in
  `operate_on_span(s)` (`field_box.hpp`, `field_tiles.hpp`); core has no SAMRAI fuser includes
- `FieldData` / `TensorFieldData::copyFuseable`: GPU_UNIFIED -> `copy_<SetEqualFuseable>`, else copy
- `Refiner::makeSchedule_` attaches the strategy (RAJA + GPU_UNIFIED fields). Schedules built
  elsewhere (synchronizers, reflux, patch ghost, regrid, MHD, diagnostics) don't, and copy as before
- test: `tests/amr/messengers/fuser_ops_strategy` (built with `withGPU AND withRaja`)
- correct (layout tests 2d 25/25) but no speedup: the schedules' cost is UM faults from host side
  ops, see Open

### Allocator / memory
- `LockedQuickPool` (`umpire_allocator.hpp`) for `PHARE::data_allocator`: umpire's QuickPool is not
  thread safe, two pools' tasks allocating at once spun forever. Recursive mutex: QuickPool calls
  its own virtuals (release from allocate, allocate/deallocate/release from coalesce)
- `soavx_chunked_appending.hpp`: reserve only when short, +5%. Refine streamed 256 particle
  chunks on top of `reserve_flat`'s .66 estimate with exact reserves: quadratic copies, and the
  umpire pool grew past 24GB at 2d 3 levels >= 80 ppc (8.6GB at 80, 767MB at 75)
- particle realloc on device: `PHARE::reserve` / `PHARE::resize` (`core/vector.hpp`)
- `res/cmake/dep/raja_umpire.cmake` links `umpire::umpire` (the plain target warned per use)

### More pools on GPU
- untiled GPU: 2 pools x 1 thread halve the per patch host work (patch per idle pool); 4 = 2 with
  1-2 patches per level

### SoAVXChunked
SoAVX sorted by cell with per cell offsets (CSR over `grow(box, 1)`), as AoSMapped is AoS + CellMap:
`arrays/particle_array_soavx_chunked.hpp`. A box is a list of contiguous particle ranges.
- `sort()`: counting sort, GPU_UNIFIED via thrust, CPU serial. Offsets valid until the next
  mutation (`sorted()`), 376 sorts per coarse step at 2d 3 levels cost 0.33 s on GPU
- exporter / selector / appender, MultiBoris, pybind layout, modules `cpp_{1_1_2,2_1_4,3_1_6}_
  SoAVXChunked_{CPU,GPU_UNIFIED}`, layout tests pass 1d/2d/3d
- SoAVX* need thrust, which only GPU builds have: `src/python3/CMakeLists.txt` skips SoAVX* lines
  without withGPU, `supported_*layouts()` (cpp_etc) and `pyphare.cpp.with_thrust()` filter them
- in the harris permutations (`tools/bench/functional/permute_harris_{2d,3d}.py`): `allocator`
  permutable, GPU_UNIFIED only if `cpp.with_gpu()`, GPU runs 1 rank (one GPU, no CUDA aware MPI)

### SoAVX (untiled, sorting based) on GPU
Every component its own column, `iCell_` per dimension; MultiBoris flat backend; exporter with
`thrust::partition`; `test_multi_updater.cpp` CPU + GPU_UNIFIED. Not a pyphare layout.

### AoSTS GPU COPY: several blocks per tile
`MultiBorisBackend::move_gpu_copy` grid y = blocks per tile, staged sums `atomicAdd` per tile:
182 -> 7.09 ns per particle.

### NVTX / scopes
- `-DPHLOP_NVTX`: every phlop scope is an NVTX range (domain "phlop"), free at run time.
  `PHARE_LOG_LEVEL=2` in `tools/cmake.sh`: level 3 is too noisy, promote what's needed to 2
- memory scopes: `-DPHARE_LOG_SCOPE_FLAGS_{1,2}=time,mem,...`, run with `PHLOP_SCOPE_MEM=peak`
- `PHARE_DEVICE_PASS` for host only code in `_PHARE_ALL_FN_` functions

## Open

### UM faults: host side ops on managed data (next)
Profile (2026-11-12, SoAVXChunked GPU, 2d 100 ppc 3 levels, `tools/bench/gpu/um_by_nvtx.py`):
784k CPU UM faults, 1.55M GPU fault pages.
- field schedules fault on the GPU, not the CPU: `average::schedules` 362k GPU pages, the
  faraday/ohm schedules 170k-255k each, ~0.55 s each. Host ops in the fills touch the fields, the
  device then faults them back. Goal: no host side op on UM data, rather than prefetching
  (`mkn::gpu::prefetch` / `prefetch_host` exist, a stopgap / measure)
- to move to the device, field side (every step):
  1. time interpolation: `FieldLinearTimeInterpolate` / `interpolateOverOverlap` (next, pointwise)
  2. refine operators: `FieldRefineOperator`, `MagneticRefinePatchStrategy`
  3. `NdArrayMask::fill2D`, `CompositeFill` (`for_each_point`)
  4. `packStreamFuseable` / `unpackStreamFuseable` (MPI), the stream buffer must be device memory
  5. `FieldBorderOpTransaction` (border sum/max) is not fuseable
- particle side (init / regrid / level ghosts): `ParticlesData<SoAVXChunked>::copy` (110k CPU
  faults), `std::vector<_, UmpireAllocator>::operator=` (column copies), the split
  (`PatternDispatcher<PurplePattern>`, `ParticlesRefineOperator`), `firstStep` / `lastStep`
- check after each: rerun nsys + `um_by_nvtx.py`, CPU faults in advance ranges -> 0

### Bugs
- **3d L2 level ghost bug (AoSTS GPU)**: `LayoutAdvanceL2Test...layout_3` fails, AoSPCTS GPU,
  AoSCMTS CPU and AoSTS CPU (debug module) pass
  - 3d: L1 ok at t=0.001, wrong at 0.00125 (first L1 step of the 2nd coarse cycle): E wrong at
    L1 x borders only, level ghosts `vx` wrong (iCells match), then L2 / L0 contaminated
  - 2d repro: seed as a parameter, AoSTS GPU fails 11/40 seeds, seed 3 in the first L1 step: L1 E
    wrong in the low-x ghost layer (~20% rel)
  - `test_multi_updater` AoSTS GPU passes all 40 seeds: the bug is on the AMR side
  - next: AoSTS CPU vs GPU at the first L1 step in 2d (seed 3), E after `fillElectricGhosts`;
    needs `cpp_2_1_4_AoSTS_CPU` and a runner monkeypatching `check_model_layout` (the old runners
    are lost)
- **AoSTS GPU 2d 3 levels segfaults** (100 ppc, 2 levels already): signal 11, null address, on a
  pool/launcher thread, first advance. Not investigated
- **4 rank harris 2d/3d at 3 levels**: SAMRAI `RefineSchedule::finishScheduleConstruction`
  internal error ("No coarser levels") for every layout incl. AoSMapped: config, not GPU

### Performance
- SoAVXChunked CPU vs AoSMapped CPU: moveIons 1.8x, `fillIonGhostParticles` 0.43 vs 0.001 s
- SoAVXChunked `initLevel` at 100 ppc still grows faster than particles (1 GB vs ~340 MB linear)
- gaps for patch ghost leavers: the push records indices leaving into the patch ghost box (atomic
  list), the exporter copies those only instead of `thrust::partition` over all; the sort stays
- Faraday / Ohm as region jobs like Ampere (low gain now), cache region jobs per level until regrid
- per call managed allocations (job lists, exporter box lists) -> persistent per level storage
- `gaps` vectors: `UmpireAllocator<std::size_t, false>`, no host zero fill
- dedupe the three transformers' GPU paths

### Build
- mkn tests on the RAJA stack: `-P samrai.ver=phare` plus `/mkn/r/llnl/samrai/phare/mkn.yaml` lib
  profile needing RAJA / camp / Umpire / cudart (headers have HAVE_RAJA/HAVE_UMPIRE)
- CPU only fusion test: a second LLNL stack without `ENABLE_CUDA` (host RAJA WorkGroup)
- SAMRAI's exported targets link the legacy `umpire` target (2 CMake warnings, fix in its build)

### Benchmark
- MultiBoris: `tools/bench/core/numerics/pusher/multi_boris.cpp` (mkn `-p bench`), results in
  `tools/bench/core/numerics/pusher/results/`; todo: `CMakeLists.txt` entry, scenario tool from a
  run's boxes / particle counts, GPU launcher threads knob

### Cleanup (debug only)
- `res/sim/gpu.txt`: `# DEBUG` + `3,1,6,AoSTS,CPU` line
- `phare_outputs/advance_hybrid/adv_layout_L2AoSTSCPU_*` (~5.6G), kept L2 GPU failure outputs
- test tolerances: `test_multi_updater` scales GPU x1e3 / tiled x1e1, layout tests don't

## How to profile (reproducible)

- harris permutations: `PHARE_SCOPE_TIMING=1 PYTHONPATH=$PWD:$PWD/build:$PWD/pyphare:<phlop>
  python3 -Ou tools/bench/functional/permute_harris_2d.py [--ppc ..] [--max_nbr_levels ..]
  [--model_layout ..] [--allocator ..] [--launches m:p:t ..]`, results `.phare_bench/<date>.tar.gz`
  (a tarball per run, timers in `.phare/timings/rank.*.txt`), `tools/bench/functional/scope_times.py`
  for per launch times. Cap memory: `systemd-run --user --scope -p MemoryMax=24G -p
  MemorySwapMax=0` (OOM kills show in `journalctl --user`), append the log (`>>`), ranks share it
- one permutation without the skip rules: `mpirun -n 1 python3 -u <script> --permutation '<json
  record>' --out-dir <dir>`, env `PHARE_THREAD_POOLS` / `PHARE_THREADS_PER_POOL`
- nsys with UM faults: `nsys profile -t cuda,nvtx,osrt --cuda-um-cpu-page-faults=true
  --cuda-um-gpu-page-faults=true -o <out> python3 ...`, `nsys export -t sqlite`, then
  `python3 tools/bench/gpu/um_by_nvtx.py <out>.sqlite [top]`: CPU faults and GPU fault pages per
  NVTX range (self / inclusive, per thread nesting) and CPU faults per faulting host function
- `nsys stats -r nvtx_gpu_proj_sum`: GPU work per NVTX range (see mkn.gpu `doc/tracing/README.md`)
- totals: `nsys stats -r cuda_gpu_kern_sum,cuda_api_sum,cuda_gpu_mem_size_sum,um_sum`
- memory blow ups: a watcher sampling the rank's VmRSS that runs `gdb -p <pid> -batch -ex
  "thread apply all bt"` past a limit (match the python rank, not mpirun); `compute-sanitizer
  --print-limit 2000` for device faults (ignore its "duplicate cub kernel" API errors)
