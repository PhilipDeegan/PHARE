# MultiBoris benchmark results, 2026-11-11

Same machine, build and run as `2026-11-10_multi_boris.md`; branch tiling_gpu at f64edda8 + uncommitted
AoSTS GPU COPY with several blocks per tile (`BlockPerLauncher` grid y, 128 threads per block,
`blocks_for(max tile particles)`, staged sums flushed with `atomicAdd` when a tile has more than one block).

- run: `PHARE_PATCHES=10 PHARE_CELLS=10 PHARE_PPC=100 ./bin/bench/phare --benchmark_min_time=0.5s` (one run, with
  launcher counters)

## Time per particle per move

| layout | alloc | COPY (push + deposit) | REF (push in place) | vs 2026-11-10 |
|---|---|---|---|---|
| AoSMapped (serial BorisPusher, master path) | CPU | 72.0 ns | 62.3 ns | = |
| AoS | CPU | 55.8 ns | 43.5 ns | = |
| AoSCMTS | CPU | 63.2 ns | 46.0 ns | = |
| AoSPCTS | CPU | 63.7 ns | 82.3 ns | = |
| SoAVX | CPU | 57.5 ns | 47.8 ns | = |
| AoSTS | GPU_UNIFIED | **7.09 ns** | 26.5 ns | COPY 168 - 182 ns -> 7.09 ns |
| AoSPCTS | GPU_UNIFIED | 92.0 ns | 83.3 ns | = |
| SoAVX | GPU_UNIFIED | 6.97 ns | 0.87 ns | = |

- AoSTS GPU COPY is now on par with SoAVX GPU COPY; the AoSTS `move_ns` counter (32 us -> 108 us) times the
  asynchronous launch only, not the kernel
- AoSTS REF vs SoAVX REF is not like for like: AoSTS REF includes the tiled tracking (part of a sort: particles
  binned into tiles), SoAVX REF does no sorting (its sort is a later step, not timed here)

## Raw output

```
2026-11-11T10:47:05+01:00
Running ./bin/bench/phare
Run on (16 X 601 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1024 KiB (x8)
  L3 Unified 98304 KiB (x1)
Load Average: 0.62, 0.40, 0.36
***WARNING*** Library was built as DEBUG. Timings may be affected.
------------------------------------------------------------------------------------------------------------------------------------------------------------
Benchmark                                                                                                  Time             CPU   Iterations UserCounters...
------------------------------------------------------------------------------------------------------------------------------------------------------------
multi_boris<dim, LayoutMode::AoSMapped, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time             72.0 ms         72.0 ms           10 items_per_second=13.8909M/s particles=1M patches=10 pools=1 t_per_particle=71.9895ns threads=1
multi_boris<dim, LayoutMode::AoSMapped, AllocatorMode::CPU, MultiBorisMode::REF>/real_time              62.3 ms         62.3 ms           11 items_per_second=16.0497M/s particles=1M patches=10 pools=1 t_per_particle=62.3064ns threads=1
multi_boris<dim, LayoutMode::AoS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time                   55.8 ms        0.089 ms           12 items_per_second=17.9282M/s particles=1M patches=10 pools=1 t_per_particle=55.778ns threads=1
multi_boris<dim, LayoutMode::AoS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                    43.5 ms        0.074 ms           16 items_per_second=23.0061M/s particles=1M patches=10 pools=1 t_per_particle=43.4668ns threads=1
multi_boris<dim, LayoutMode::AoSCMTS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time               63.2 ms        0.158 ms           11 items_per_second=15.8278M/s particles=1M patches=10 pools=1 t_per_particle=63.1801ns threads=1
multi_boris<dim, LayoutMode::AoSCMTS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                46.0 ms        0.114 ms           15 items_per_second=21.7489M/s particles=1M patches=10 pools=1 t_per_particle=45.9792ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time               63.7 ms        0.191 ms           11 items_per_second=15.7025M/s particles=1M patches=10 pools=1 t_per_particle=63.6841ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                82.3 ms        0.224 ms            9 items_per_second=12.1475M/s particles=1M patches=10 pools=1 t_per_particle=82.3216ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time                 57.5 ms        0.117 ms           12 items_per_second=17.3996M/s particles=1M patches=10 pools=1 t_per_particle=57.4726ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                  47.8 ms        0.103 ms           15 items_per_second=20.9244M/s particles=1M patches=10 pools=1 t_per_particle=47.7912ns threads=1
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time         7.09 ms         7.03 ms           85 items_per_second=141.13M/s move_ns=108.208k particles=1M patches=10 pools=1 sync_ns=372.052 t_per_particle=7.08564ns threads=1
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time          26.5 ms         21.8 ms           27 items_per_second=37.7979M/s move_ns=16.7143M particles=1M patches=10 pools=1 sync_ns=8.78268M t_per_particle=26.4565ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time       92.0 ms         91.6 ms            8 items_per_second=10.8703M/s move_ns=54.7282M particles=1M patches=10 pools=1 sync_ns=356.2 t_per_particle=91.9939ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time        83.3 ms         79.1 ms           12 items_per_second=12.0075M/s move_ns=35.2507M particles=1M patches=10 pools=1 sync_ns=46.3183M t_per_particle=83.281ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time         6.97 ms         6.74 ms           85 items_per_second=143.566M/s move_ns=3.85294M particles=1M patches=10 pools=1 sync_ns=362.174 t_per_particle=6.96544ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time         0.873 ms        0.832 ms          663 items_per_second=1.14491G/s move_ns=424.14k particles=1M patches=10 pools=1 sync_ns=1.34933k t_per_particle=873.429ps threads=1
```
