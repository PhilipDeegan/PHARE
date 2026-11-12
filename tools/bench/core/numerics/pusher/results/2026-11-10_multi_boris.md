# MultiBoris benchmark results, 2026-11-10

`tools/bench/core/numerics/pusher/multi_boris.cpp`, see doc/gpu.md for what is timed.

- machine: AMD Ryzen 7 9800X3D (8 cores / 16 threads), NVIDIA GeForce RTX 5080 (16 GB, driver 620.70.06, sm_120)
- branch tiling_gpu at e8eaa27b + uncommitted changes (SoAVX GPU path, transformers, ...)
- build: `mkn clean build -p bench -M tools/bench/core/numerics/pusher/multi_boris.cpp -w mkn.gpu -x res/mkn/clang_cuda -O 3 -a -DNDEBUG`
  (google benchmark itself is a debug build, its warning below; doesn't affect what is timed)
- run: `PHARE_PATCHES=10 PHARE_CELLS=10 PHARE_PPC=100 ./bin/bench/phare --benchmark_min_time=0.5s`
  synthetic level: 10 patches x 10^3 cells x 100 ppc = 1M particles (+ a 1 cell level ghost layer each), 3d,
  v = 0, defaults otherwise: CPU 1 pool x 1 thread, GPU 1 ThreadedStreamLauncher host thread + a stream per patch

## Time per particle per move

| layout | alloc | COPY (push + deposit) | REF (push in place) |
|---|---|---|---|
| AoSMapped (serial BorisPusher, master path) | CPU | 72.4 ns | 62.3 ns |
| AoS | CPU | 56.0 ns | 43.7 ns |
| AoSCMTS | CPU | 63.4 ns | 46.4 ns |
| AoSPCTS | CPU | 63.8 ns | 81.8 ns |
| SoAVX | CPU | 57.4 ns | 47.7 ns |
| AoSTS | GPU_UNIFIED | 168 - 182 ns | 25.3 - 26.1 ns |
| AoSPCTS | GPU_UNIFIED | 92.0 ns | 82.1 - 83.9 ns |
| SoAVX | GPU_UNIFIED | 6.6 - 6.7 ns | 0.78 - 0.90 ns |

GPU ranges: the two runs below (the second with the launcher time counters).

## GPU launcher times, per patch per iteration (second run)

Intervals from a function's start on its patch to its end, waits on the patch's stream included, patches
overlap: per patch, not additive. "move": push (REF) or COPY's launch, "sync": REF stream wait + sync_ref.

| layout | mode | move | sync | reading |
|---|---|---|---|---|
| AoSTS | COPY | 32 us | 0.6 us | async launch only timed: the time is the tile deposit kernel |
| AoSTS | REF | 16.5 ms | 8.6 ms | push + tiled tracking |
| AoSPCTS | COPY | 54.7 ms | 0.4 us | no GPU copy kernel: COPY on the host over unified memory |
| AoSPCTS | REF | 34.7 ms | 44.9 ms | sync_ref (tracking) the larger part |
| SoAVX | COPY | 3.65 ms | 0.4 us | push + deposit kernel per patch (launch waits its stream) |
| SoAVX | REF | 0.44 ms | 1.3 us | push kernel |

## Raw output

### CPU + GPU
```
2026-11-10T19:25:36+01:00
Running ./bin/bench/phare
Run on (16 X 4353.62 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1024 KiB (x8)
  L3 Unified 98304 KiB (x1)
Load Average: 1.64, 2.25, 2.20
***WARNING*** Library was built as DEBUG. Timings may be affected.
------------------------------------------------------------------------------------------------------------------------------------------------------------
Benchmark                                                                                                  Time             CPU   Iterations UserCounters...
------------------------------------------------------------------------------------------------------------------------------------------------------------
multi_boris<dim, LayoutMode::AoSMapped, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time             72.4 ms         72.4 ms           10 items_per_second=13.8197M/s particles=1M patches=10 pools=1 t_per_particle=72.3606ns threads=1
multi_boris<dim, LayoutMode::AoSMapped, AllocatorMode::CPU, MultiBorisMode::REF>/real_time              62.3 ms         62.3 ms           11 items_per_second=16.0457M/s particles=1M patches=10 pools=1 t_per_particle=62.3219ns threads=1
multi_boris<dim, LayoutMode::AoS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time                   56.0 ms        0.098 ms           12 items_per_second=17.8557M/s particles=1M patches=10 pools=1 t_per_particle=56.0045ns threads=1
multi_boris<dim, LayoutMode::AoS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                    43.7 ms        0.113 ms           16 items_per_second=22.9024M/s particles=1M patches=10 pools=1 t_per_particle=43.6635ns threads=1
multi_boris<dim, LayoutMode::AoSCMTS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time               63.4 ms        0.167 ms           11 items_per_second=15.7851M/s particles=1M patches=10 pools=1 t_per_particle=63.351ns threads=1
multi_boris<dim, LayoutMode::AoSCMTS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                46.4 ms        0.137 ms           15 items_per_second=21.5442M/s particles=1M patches=10 pools=1 t_per_particle=46.4163ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time               63.8 ms        0.234 ms           11 items_per_second=15.6637M/s particles=1M patches=10 pools=1 t_per_particle=63.8418ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                81.8 ms        0.230 ms            9 items_per_second=12.2289M/s particles=1M patches=10 pools=1 t_per_particle=81.7738ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::CPU, MultiBorisMode::COPY>/real_time                 57.4 ms        0.126 ms           12 items_per_second=17.4328M/s particles=1M patches=10 pools=1 t_per_particle=57.3632ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::CPU, MultiBorisMode::REF>/real_time                  47.7 ms        0.120 ms           15 items_per_second=20.949M/s particles=1M patches=10 pools=1 t_per_particle=47.7349ns threads=1
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time          168 ms          168 ms            4 items_per_second=5.94M/s particles=1M patches=10 pools=1 t_per_particle=168.35ns threads=1
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time          25.3 ms         20.3 ms           27 items_per_second=39.567M/s particles=1M patches=10 pools=1 t_per_particle=25.2736ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time       91.9 ms         91.4 ms            8 items_per_second=10.8803M/s particles=1M patches=10 pools=1 t_per_particle=91.9093ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time        83.9 ms         79.1 ms           13 items_per_second=11.9124M/s particles=1M patches=10 pools=1 t_per_particle=83.9463ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time         6.73 ms         6.47 ms           93 items_per_second=148.607M/s particles=1M patches=10 pools=1 t_per_particle=6.72917ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time         0.784 ms        0.748 ms          696 items_per_second=1.27511G/s particles=1M patches=10 pools=1 t_per_particle=784.246ps threads=1
```

### GPU, with launcher counters
```
2026-11-10T19:34:40+01:00
Running ./bin/bench/phare
Run on (16 X 4170.41 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1024 KiB (x8)
  L3 Unified 98304 KiB (x1)
Load Average: 3.04, 2.61, 2.34
***WARNING*** Library was built as DEBUG. Timings may be affected.
------------------------------------------------------------------------------------------------------------------------------------------------------------
Benchmark                                                                                                  Time             CPU   Iterations UserCounters...
------------------------------------------------------------------------------------------------------------------------------------------------------------
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time          182 ms          182 ms            4 items_per_second=5.50824M/s move_ns=31.6216k particles=1M patches=10 pools=1 sync_ns=584.65 t_per_particle=181.546ns threads=1
multi_boris<dim, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time          26.1 ms         21.6 ms           26 items_per_second=38.3221M/s move_ns=16.4596M particles=1M patches=10 pools=1 sync_ns=8.63021M t_per_particle=26.0946ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time       92.0 ms         91.5 ms            8 items_per_second=10.8736M/s move_ns=54.6944M particles=1M patches=10 pools=1 sync_ns=358.637 t_per_particle=91.9658ns threads=1
multi_boris<dim, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time        82.1 ms         75.5 ms           12 items_per_second=12.1758M/s move_ns=34.7312M particles=1M patches=10 pools=1 sync_ns=44.9418M t_per_particle=82.13ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::COPY>/real_time         6.64 ms         6.41 ms           91 items_per_second=150.508M/s move_ns=3.6547M particles=1M patches=10 pools=1 sync_ns=355.805 t_per_particle=6.64416ns threads=1
multi_boris<dim, LayoutMode::SoAVX, AllocatorMode::GPU_UNIFIED, MultiBorisMode::REF>/real_time         0.900 ms        0.863 ms          627 items_per_second=1.11171G/s move_ns=440.98k particles=1M patches=10 pools=1 sync_ns=1.30551k t_per_particle=899.511ps threads=1
```
