# UM page faults attributed to phlop NVTX ranges (nsys export -t sqlite)
#  usage: python3 um_by_nvtx.py <file.sqlite> [top]
#  CPU fault rows have no thread: attribution is by time only, inclusive counts per range via
#  searchsorted, self = inclusive - children (per thread nesting)
import collections
import sqlite3
import sys

import numpy as np

db = sqlite3.connect(sys.argv[1])
top = int(sys.argv[2]) if len(sys.argv) > 2 else 25
strings = dict(db.execute("select id, value from StringIds"))


def has(table):
    return db.execute(
        "select count(*) from sqlite_master where type='table' and name=?", (table,)
    ).fetchone()[0]


rows = db.execute(
    "select start, end, globalTid, coalesce(text, ''), textId from NVTX_EVENTS "
    "where eventType in (59, 60) and end is not null"
).fetchall()
ranges = []
for start, end, tid, text, text_id in rows:
    name = text or strings.get(text_id, "?")
    ranges.append((start, end, tid, name))
ranges.sort(key=lambda r: (r[2], r[0], -r[1]))
print(f"{len(ranges)} NVTX ranges")

cpu = np.sort(
    np.array([r[0] for r in db.execute("select start from CUDA_UM_CPU_PAGE_FAULT_EVENTS")])
    if has("CUDA_UM_CPU_PAGE_FAULT_EVENTS")
    else np.array([])
)
gpu_rows = (
    db.execute("select start, numberOfPageFaults from CUDA_UM_GPU_PAGE_FAULT_EVENTS").fetchall()
    if has("CUDA_UM_GPU_PAGE_FAULT_EVENTS")
    else []
)
gpu_rows.sort()
gpu_t = np.array([r[0] for r in gpu_rows])
gpu_cum = np.concatenate([[0], np.cumsum([r[1] for r in gpu_rows])])
print(f"{len(cpu)} CPU UM faults, {gpu_cum[-1]} GPU fault pages")


def count(t, a, b):
    return np.searchsorted(t, b) - np.searchsorted(t, a)


def gpu_count(a, b):
    return gpu_cum[np.searchsorted(gpu_t, b)] - gpu_cum[np.searchsorted(gpu_t, a)]


incl = [(count(cpu, s, e), gpu_count(s, e)) for s, e, _, _ in ranges]

# per thread nesting: parent = innermost open range on the same thread
children = collections.defaultdict(lambda: [0, 0])
stack = []
for i, (s, e, tid, _) in enumerate(ranges):
    while stack and (ranges[stack[-1]][2] != tid or ranges[stack[-1]][1] <= s):
        stack.pop()
    if stack:
        p = stack[-1]
        children[p][0] += incl[i][0]
        children[p][1] += incl[i][1]
    stack.append(i)

agg = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])  # calls, cpu incl, cpu self, gpu incl, gpu self, ns
for i, (s, e, _, name) in enumerate(ranges):
    a = agg[name]
    a[0] += 1
    a[1] += incl[i][0]
    a[2] += incl[i][0] - children[i][0]
    a[3] += incl[i][1]
    a[4] += incl[i][1] - children[i][1]
    a[5] += e - s

print(f"\n{'range (self = not in a child range)':58}{'calls':>6}{'cpu self':>10}{'cpu incl':>10}{'gpu self':>10}{'time s':>8}")
for name, a in sorted(agg.items(), key=lambda kv: -kv[1][2])[:top]:
    print(f"{name[:58]:58}{a[0]:6}{a[2]:10}{a[1]:10}{a[4]:10}{a[5] / 1e9:8.3f}")

if has("CUDA_UM_CPU_PAGE_FAULT_EVENTS"):
    print("\nCPU UM faults by faulting host function:")
    for fn, c in db.execute(
        "select CpuInstruction, count(*) c from CUDA_UM_CPU_PAGE_FAULT_EVENTS "
        "group by CpuInstruction order by c desc limit 15"
    ):
        print(f"{c:9}  {strings.get(fn, fn)}"[:200])
