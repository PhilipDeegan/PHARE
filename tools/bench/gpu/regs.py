#!/usr/bin/env python3
"""
Per-kernel GPU resource usage of built binaries, and the largest block size
each kernel can launch with. mkn.gpu fits guessed launches to this limit at
runtime; Fixed launches (caller chosen block) must stay under it themselves.

    tools/bench/gpu/regs.py                                   # bin/test_core/phare
    tools/bench/gpu/regs.py -t 512                            # kernels that can't launch at 512
    tools/bench/gpu/regs.py dbg=bin/a/phare rel=bin/b/phare   # compare build configurations

A kernel using R registers per thread gets them per warp, rounded up to the
allocation unit; a block of T threads needs ceil(T/32) such warps out of the
registers-per-block budget. Launching more threads than fit fails with
"too many resources requested for launch".
"""

import argparse
import re
import shutil
import subprocess
import sys

WARP = 32

FN_RE = re.compile(r"Function (\S+):")
RES_RE = re.compile(r"REG:(\d+) STACK:(\d+) SHARED:(\d+) LOCAL:(\d+)")
ARCH_RE = re.compile(r"arch = (sm_\d+)")


def cuobjdump():
    exe = shutil.which("cuobjdump") or "/usr/local/cuda/bin/cuobjdump"
    return exe


def demangle(names):
    if not names:
        return {}
    # binutils c++filt leaves C++20 mangling (e.g. auto NTTPs) untouched, llvm's handles it
    exe = shutil.which("llvm-cxxfilt") or "c++filt"
    out = subprocess.run(
        [exe], input="\n".join(names), capture_output=True, text=True
    ).stdout.splitlines()
    return dict(zip(names, out))


def shorten(name, depth):
    """drop template arguments nested deeper than depth"""
    out, level = [], 0
    for c in name:
        if c == "<":
            level += 1
            if level == depth + 1:
                out.append("<…>")
            if level > depth:
                continue
        elif c == ">":
            level -= 1
            if level >= depth:
                continue
        elif level > depth:
            continue
        out.append(c)
    return "".join(out)


def read_kernels(binary):
    res = subprocess.run(
        [cuobjdump(), "--dump-resource-usage", binary], capture_output=True, text=True
    )
    if res.returncode != 0:
        sys.exit(f"cuobjdump failed on {binary}: {res.stderr.strip()}")

    kernels, arch, fn = [], "?", None
    for line in res.stdout.splitlines():
        if m := ARCH_RE.search(line):
            arch = m.group(1)
        if m := FN_RE.search(line):
            fn = m.group(1)
        if fn and (m := RES_RE.search(line)):
            reg, stack, shared, local = map(int, m.groups())
            kernels.append(dict(arch=arch, mangled=fn, reg=reg, stack=stack,
                                shared=shared, local=local))  # fmt: skip
            fn = None

    names = demangle([k["mangled"] for k in kernels])
    for k in kernels:
        k["name"] = names.get(k["mangled"], k["mangled"])
    return kernels


def max_threads(reg, regs_per_block, alloc_unit, max_block):
    if reg == 0:
        return max_block
    per_warp = -(-reg * WARP // alloc_unit) * alloc_unit  # ceil to allocation unit
    return min(max_block, (regs_per_block // per_warp) * WARP)


def pow2_floor(n):
    return 1 << (n.bit_length() - 1) if n > 0 else 0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("binaries", nargs="*", default=["bin/test_core/phare"],
                   help="path or label=path (default bin/test_core/phare)")  # fmt: skip
    p.add_argument("-t", "--threads", type=int,
                   help="only kernels that can't launch with this many threads per block")  # fmt: skip
    p.add_argument("-n", "--top", type=int, default=20, help="rows per binary (0 = all)")
    p.add_argument("-d", "--depth", type=int, default=1, help="template depth shown")
    p.add_argument("--filter", default="", help="regex on the demangled kernel name")
    p.add_argument("--regs-per-block", type=int, default=65536)
    p.add_argument("--alloc-unit", type=int, default=256,
                   help="registers allocated per warp in multiples of this")  # fmt: skip
    p.add_argument("--max-block", type=int, default=1024)
    args = p.parse_args()

    limit = lambda reg: max_threads(  # noqa: E731
        reg, args.regs_per_block, args.alloc_unit, args.max_block
    )
    flt = re.compile(args.filter) if args.filter else None

    summary = {}
    for spec in args.binaries:
        label, _, path = spec.rpartition("=")
        label = label or path
        kernels = sorted(read_kernels(path), key=lambda k: -k["reg"])
        if flt:
            kernels = [k for k in kernels if flt.search(k["name"])]

        worst = min((limit(k["reg"]) for k in kernels), default=args.max_block)
        summary[label] = (len(kernels), max((k["reg"] for k in kernels), default=0), worst)

        rows = kernels
        if args.threads:
            rows = [k for k in kernels if limit(k["reg"]) < args.threads]
        if args.top:
            rows = rows[: args.top]

        print(f"\n== {label} ({path}): {len(kernels)} kernels")
        print(f"{'REG':>4} {'STACK':>6} {'SHARED':>6} {'MAXTHR':>6}  arch     kernel")
        for k in rows:
            print(f"{k['reg']:>4} {k['stack']:>6} {k['shared']:>6} {limit(k['reg']):>6}  "
                  f"{k['arch']:<8} {shorten(k['name'], args.depth)}")  # fmt: skip

    print("\n== summary: largest block size every kernel can launch with")
    print(f"{'kernels':>8} {'maxREG':>6} {'maxthr':>6} {'pow2':>5}  config")
    for label, (n, reg, worst) in summary.items():
        print(f"{n:>8} {reg:>6} {worst:>6} {pow2_floor(worst):>5}  {label}")



if __name__ == "__main__":
    main()
