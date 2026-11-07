import dataclasses
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path

FILE_DIR = Path(__file__).resolve().parent
PHARE_BUILD_DIR = Path(os.environ["PHARE_BUILD_DIR"])
BUILD_DIR = PHARE_BUILD_DIR / "phare_configurator"  # see cmake.sh
PROBE_DIR = BUILD_DIR / "probe"
GPU_DIR = BUILD_DIR / "gpu"
LOCAL_CMAKE = BUILD_DIR / "local.cmake"  # included by res/cmake/def.cmake
ROOT_DIR = FILE_DIR.parent.parent
DOT_PHARE_DIR = ROOT_DIR / ".phare"
FILES = list(BUILD_DIR.glob("PHARE_*"))
DEF_DIR = PHARE_BUILD_DIR / "src" / "core" / "def"
GENERATED_CONFIGS = dict(
    mpi=DEF_DIR / "_gen_mpi.hpp",
    system=DEF_DIR / "_gen_sys.hpp",
)


@dataclasses.dataclass
class SystemSettings:
    cmake_binary: str
    cmake_version: str
    cxx_compiler: str
    cxx_compiler_version: str
    hdf5_version: str
    hdf5_is_parallel: str
    mpi_version: str
    python_binary: str
    python_version: str
    uname: str
    git_hash: str
    git_message: str
    gpu_vendor: str
    gpu_arch: str
    gpu_devices: str


SYSTEM_CPP_ = """
#include <string_view>
#include <unordered_map>

namespace PHARE {{
struct SystemConfig {{

    constexpr static std::string_view CMAKE_BINARY = R"({})";
    constexpr static std::string_view CMAKE_VERSION = R"({})";
    constexpr static std::string_view CXX_COMPILER = R"({})";
    constexpr static std::string_view CXX_COMPILER_VERSION = R"({})";
    constexpr static std::string_view HDF5_VERSION = R"({})";
    constexpr static std::string_view HDF5_IS_PARALLEL = R"({})";
    constexpr static std::string_view _MPI_VERSION_ = R"({})";
    constexpr static std::string_view PYTHON_BINARY = R"({})";
    constexpr static std::string_view PYTHON_VERSION = R"({})";
    constexpr static std::string_view UNAME = R"({})";
    constexpr static std::string_view GIT_HASH = R"({})";
    constexpr static std::string_view GIT_MESSAGE = R"({})";
    constexpr static std::string_view GPU_VENDOR = R"({})";
    constexpr static std::string_view GPU_ARCH = R"({})";
    constexpr static std::string_view GPU_DEVICES = R"({})";

}};

inline std::unordered_map<std::string, std::string> build_config(){{
  return {{
      {{"CMAKE_BINARY", std::string{{SystemConfig::CMAKE_BINARY}}}},
      {{"CMAKE_VERSION", std::string{{SystemConfig::CMAKE_VERSION}}}},
      {{"CXX_COMPILER", std::string{{SystemConfig::CXX_COMPILER}}}},
      {{"CXX_COMPILER_VERSION", std::string{{SystemConfig::CXX_COMPILER_VERSION}}}},
      {{"HDF5_VERSION", std::string{{SystemConfig::HDF5_VERSION}}}},
      {{"HDF5_IS_PARALLEL", std::string{{SystemConfig::HDF5_IS_PARALLEL}}}},
      {{"MPI_VERSION", std::string{{SystemConfig::_MPI_VERSION_}}}},
      {{"PYTHON_BINARY", std::string{{SystemConfig::PYTHON_BINARY}}}},
      {{"PYTHON_VERSION", std::string{{SystemConfig::PYTHON_VERSION}}}},
      {{"UNAME", std::string{{SystemConfig::UNAME}}}},
      {{"GIT_HASH", std::string{{SystemConfig::GIT_HASH}}}},
      {{"GIT_MESSAGE", std::string{{SystemConfig::GIT_MESSAGE}}}},
      {{"GPU_VENDOR", std::string{{SystemConfig::GPU_VENDOR}}}},
      {{"GPU_ARCH", std::string{{SystemConfig::GPU_ARCH}}}},
      {{"GPU_DEVICES", std::string{{SystemConfig::GPU_DEVICES}}}}
  }};
}}

}}

"""


def subprocess_run(cmd, on_error="error message"):
    try:
        proc = subprocess.run(cmd.split(" "), check=True, capture_output=True)
        return (proc.stdout + proc.stderr).decode().strip()
    except (OSError, subprocess.SubprocessError):
        return on_error


def file_string_or(filename, fail=None):
    if not fail:
        fail = f"couldn't open file {filename}"
    filepath = BUILD_DIR / filename
    if not filepath.exists():
        return fail
    with open(filepath) as f:
        return f.readline().strip()


def cmake_cache_line_to_key_value_pair(line):
    bits = line.split("=")
    return bits[0], "=".join(bits[1:])


def parse_cmake_cache_file():
    cmaka_cache_file = PROBE_DIR / "CMakeCache.txt"
    dic = {}
    with open(cmaka_cache_file) as f:
        for line in f.readlines():
            line = line.strip()
            if any([line.startswith("#"), "=" not in line]):
                continue
            key, val = cmake_cache_line_to_key_value_pair(line)
            dic[key] = val
    return dic


def cmake_version_from(cmake_cache: dict) -> str:
    if all(
        [
            "CMAKE_CACHE_MAJOR_VERSION:INTERNAL" in cmake_cache,
            "CMAKE_CACHE_MINOR_VERSION:INTERNAL" in cmake_cache,
            "CMAKE_CACHE_PATCH_VERSION:INTERNAL" in cmake_cache,
        ]
    ):
        major = cmake_cache["CMAKE_CACHE_MAJOR_VERSION:INTERNAL"]
        minor = cmake_cache["CMAKE_CACHE_MINOR_VERSION:INTERNAL"]
        patch = cmake_cache["CMAKE_CACHE_PATCH_VERSION:INTERNAL"]
        return f"{major}.{minor}.{patch}"

    return "Failed to find version in cache file"


def local_file_format(filename):
    return filename[6:-4]


def create_file(f):
    with open(f, "a") as file:
        file.truncate(0)


def ifndef_define(var, val, buf, comment=None):
    buf += f"#ifndef {var} \n"
    if comment:
        buf += f"// {comment} \n"
    buf += f"#define {var} {val}\n"
    buf += f"#endif /* {var} */ \n\n"
    return buf


def deduce_mpi_type_from(version):
    if any(s in version for s in ["Open MPI", "OpenMPI"]):
        return "OMPI"
    if any(s in version for s in ["MPICH"]):
        return "MPICH"
    return "unknown"


def config_mpi_version(txtfile, h_file):
    with open(txtfile) as f:
        version = f.readline()
    buf = "\n"
    mpi_type = deduce_mpi_type_from(version)
    if mpi_type == "OMPI":
        buf = ifndef_define(
            "OMPI_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    elif mpi_type == "MPICH":
        buf = ifndef_define(
            "MPICH_SKIP_MPICXX", 1, buf, "avoids default including mpicxx"
        )
    with open(h_file, "w") as f:
        f.write(buf)
    return version


def try_compiler_dash_v(compiler):
    return subprocess_run(f"{compiler} -v", "compiler does not support '-v' argument")


def get_python_version(py3_binary):
    return subprocess_run(
        f"{py3_binary} -V", "python3 interpreter does not support '-V' argument"
    )


def gen_system_file(gpu):
    out_file = GENERATED_CONFIGS["system"]
    cmake_cache: dict = parse_cmake_cache_file()

    settings = SystemSettings(
        cmake_binary=os.environ["CMAKE_COMMAND"],
        cmake_version=cmake_version_from(cmake_cache),
        cxx_compiler=os.environ["CMAKE_CXX_COMPILER"],
        cxx_compiler_version=try_compiler_dash_v(os.environ["CMAKE_CXX_COMPILER"]),
        hdf5_version=file_string_or("PHARE_HDF5_version.txt", "HDF5 version failed"),
        hdf5_is_parallel=file_string_or(
            "PHARE_HDF5_is_parallel.txt", "HDF5 is parallel check failed"
        ),
        mpi_version=file_string_or(
            "PHARE_MPI_Get_library_version.txt", "MPI version failed"
        ),
        python_binary=os.environ["PYTHON_EXECUTABLE"],
        python_version=get_python_version(os.environ["PYTHON_EXECUTABLE"]),
        uname=subprocess_run("uname -a"),
        git_hash=subprocess_run("git log -1 --pretty=format:%H", "git log failed"),
        git_message=subprocess_run("git log -1", "git log failed"),
        gpu_vendor=gpu.vendor,
        gpu_arch=gpu.arch,
        gpu_devices="\n".join(gpu.devices),
    )
    with open(DOT_PHARE_DIR / "build_config.json", "w") as f:
        json.dump(dataclasses.asdict(settings), f)

    with open(out_file, "w") as f:
        f.write(
            SYSTEM_CPP_.format(
                *[getattr(settings, v.name) for v in dataclasses.fields(settings)]
            )
        )


def config_mpi():
    h_file = GENERATED_CONFIGS["mpi"]
    create_file(h_file)
    operators = dict(MPI_Get_library_version=config_mpi_version)
    results = {}
    for f in FILES:
        local_name = local_file_format(f.name)
        if local_name in operators:
            results[local_name] = operators[local_name](f, h_file)
    return results


def write_local_cmake_file(mpi_results, gpu):
    with open(LOCAL_CMAKE, "w") as file:
        file.write("""cmake_minimum_required (VERSION 3.20.1)
project(configured_phare)
message("")
message("!!PHARE CONFIGURATED!!")
message("")

""")

        mpi_type = deduce_mpi_type_from(mpi_results["MPI_Get_library_version"])
        if mpi_type == "OMPI":
            # work around for https://github.com/open-mpi/ompi/issues/10761#issuecomment-1236909802
            file.write(
                """set (PHARE_MPIRUN_POSTFIX "${PHARE_MPIRUN_POSTFIX};--bind-to;none")
                """
            )

        file.write(gpu.cmake())


# --- GPU ---
# PHARE_GPU: auto (default, CUDA then HIP) | off | cuda | hip
# PHARE_GPU_ARCH: e.g. sm_90 / gfx90a, skips the device query (build nodes without a GPU),
#  the toolchain check then only compiles, as there may be no device to run on
# found: vendor + arch from a host only runtime query (gpu_probe.cpp), then gpu_check.cpp is
#  built with the flags PHARE would get and run, only then are they written to local.cmake


@dataclasses.dataclass
class GpuConfig:
    vendor: str = ""
    arch: str = ""
    devices: list = dataclasses.field(default_factory=list)
    cxx_flags: list = dataclasses.field(default_factory=list)
    link_flags: list = dataclasses.field(default_factory=list)
    status: str = "not configured"

    def cmake(self):
        def cm_list(flags):
            # option pairs as one "SHELL:" item, else CMake's option de-duplication can drop
            #  a repeated "-x" / "-Xclang"
            items, it = [], iter(flags)
            for flag in it:
                items.append(
                    f"SHELL:{flag} {next(it)}" if flag in ("-x", "-Xclang") else flag
                )
            return ";".join(items)

        return f"""
# GPU, see tools/config/config.py: {self.status}
set (PHARE_GPU_VENDOR "{self.vendor}")
set (PHARE_GPU_ARCH "{self.arch}")
set (PHARE_GPU_CXX_FLAGS "{cm_list(self.cxx_flags)}")
set (PHARE_GPU_LINK_FLAGS "{cm_list(self.link_flags)}")
set (PHARE_GPU_STATUS "{self.status}")
message(STATUS "PHARE configurator GPU: {self.status}")
"""


def run_cmd(cmd, cwd=None):
    """list cmd, returns (ok, output)"""
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, capture_output=True, timeout=300, check=False
        )
        out = (proc.stdout + proc.stderr).decode(errors="replace").strip()
        return proc.returncode == 0, out
    except (OSError, subprocess.SubprocessError) as e:
        return False, str(e)


def first_error(out):
    lines = out.splitlines()
    return next((line for line in lines if "error" in line), lines[-1] if lines else "")


def first_existing(*paths):
    for path in paths:
        if path and Path(path).exists():
            return Path(path)
    return None


def cuda_toolkit():
    nvcc = shutil.which("nvcc")
    root = first_existing(
        os.environ.get("CUDA_PATH"),
        os.environ.get("CUDA_HOME"),
        Path(nvcc).resolve().parent.parent if nvcc else None,
        "/usr/local/cuda",
    )
    if not root:
        return None
    target = root / "targets" / f"{platform.machine()}-linux"  # toolkit cross layout
    inc = first_existing(root / "include", target / "include")
    lib = first_existing(root / "lib64", target / "lib", root / "lib")
    if not inc or not lib or not (inc / "cuda_runtime_api.h").exists():
        return None
    return dict(root=root, inc=inc, lib=lib, cccl=first_existing(inc / "cccl"))


def rocm_toolkit():
    hipconfig = shutil.which("hipconfig")
    ok, out = run_cmd([hipconfig, "--rocmpath"]) if hipconfig else (False, "")
    root = first_existing(
        os.environ.get("ROCM_PATH"),
        os.environ.get("HIP_PATH"),
        out.strip() if ok else None,
        "/opt/rocm",
    )
    if not root or not (root / "include/hip/hip_runtime_api.h").exists():
        return None
    return dict(root=root, inc=root / "include", lib=root / "lib")


def probe_devices(cxx, vendor, tk):
    """host only query: plain C++ compiler + runtime library"""
    exe = GPU_DIR / f"gpu_probe_{vendor}"
    cmd = [cxx, "-O1", str(FILE_DIR / "gpu_probe.cpp"), "-o", str(exe)]
    cmd += [f"-I{tk['inc']}", f"-L{tk['lib']}", f"-Wl,-rpath,{tk['lib']}"]
    if vendor == "hip":
        cmd += ["-DPHARE_PROBE_HIP=1", "-D__HIP_PLATFORM_AMD__", "-lamdhip64"]
    else:
        cmd += ["-lcudart"]
    ok, out = run_cmd(cmd)
    if not ok:
        return None, f"{vendor} probe build failed: {first_error(out)}"
    ok, out = run_cmd([str(exe)])
    if not ok or not out:
        return None, f"no usable {vendor} device"
    devices = [
        dict(kv.split("=", 1) for kv in line.split("\t") if "=" in kv)
        for line in out.splitlines()
    ]
    return devices, ""


def gpu_flags(cxx, vendor, arch, tk):
    """the flags PHARE gets, as res/mkn/clang_cuda.yaml / hip.yaml (minus std/fPIC)
    clang style (-x cuda / -x hip, no nvcc): whether cxx takes them is left to check_toolchain"""
    if vendor == "cuda":
        cxx_flags = [
            "-x", "cuda", f"--cuda-gpu-arch={arch}", f"--cuda-path={tk['root']}",
            "-Wno-unknown-cuda-version", "-Xclang", "-fcuda-allow-variadic-functions",
            "-DMKN_GPU_CUDA=1", "-DMKN_GPU_ROCM=0", "-DPHARE_HAVE_MKN_GPU=1",
        ]  # fmt: skip
        if tk["cccl"]:
            cxx_flags += [f"-I{tk['cccl']}"]
        link_flags = [f"-L{tk['lib']}", f"-Wl,-rpath,{tk['lib']}", "-lcudart"]
    else:
        cxx_flags = [
            "-x", "hip", f"--offload-arch={arch}",
            "-DMKN_GPU_ROCM=1", "-DMKN_GPU_CUDA=0", "-DPHARE_HAVE_MKN_GPU=1",
        ]  # fmt: skip
        link_flags = [f"-L{tk['lib']}", f"-Wl,-rpath,{tk['lib']}", "-lamdhip64"]
    return cxx_flags, link_flags


def check_toolchain(cxx, vendor, cxx_flags, link_flags, run):
    exe = GPU_DIR / f"gpu_check_{vendor}"
    hip = ["-DPHARE_PROBE_HIP=1"] if vendor == "hip" else []
    cmd = [cxx, *cxx_flags, *hip, str(FILE_DIR / "gpu_check.cpp"), "-o", str(exe)]
    ok, out = run_cmd(cmd + link_flags)
    if not ok:
        return f"{cxx} cannot build for it: {first_error(out)}"
    if run:
        ok, out = run_cmd([str(exe)])
        if not ok:
            return f"test kernel failed to run {out}"
    return ""


def config_gpu():
    mode = os.environ.get("PHARE_GPU", "auto").lower() or "auto"
    forced_arch = os.environ.get("PHARE_GPU_ARCH", "")
    cxx = os.environ["CMAKE_CXX_COMPILER"]

    if mode in ("off", "0", "false", "no"):
        return GpuConfig(status="disabled (PHARE_GPU=off)")

    vendors = ["cuda", "hip"] if mode == "auto" else [mode]
    toolkits = dict(cuda=cuda_toolkit, hip=rocm_toolkit)
    reasons = []
    for vendor in vendors:
        if vendor not in toolkits:
            return GpuConfig(status=f"unknown PHARE_GPU={mode}")
        tk = toolkits[vendor]()
        if not tk:
            reasons.append(f"no {vendor} toolkit")
            continue

        devices = []
        if forced_arch:
            arch = forced_arch
        else:
            devices, why = probe_devices(cxx, vendor, tk)
            if not devices:
                reasons.append(why)
                continue
            arch = devices[0].get("arch", "")
            if any(d.get("arch") != arch for d in devices):
                reasons.append(f"mixed {vendor} archs, using device 0 ({arch})")

        cxx_flags, link_flags = gpu_flags(cxx, vendor, arch, tk)
        why = check_toolchain(cxx, vendor, cxx_flags, link_flags, run=bool(devices))
        if why:
            reasons.append(f"{vendor} {arch}: {why}")
            continue

        described = ["\t".join(f"{k}={v}" for k, v in d.items()) for d in devices] or [
            f"arch={arch} (PHARE_GPU_ARCH, not run)"
        ]
        status = f"{vendor} {arch}" + (
            f" ({devices[0].get('name')})" if devices else ""
        )
        if reasons:
            status += " [" + "; ".join(reasons) + "]"
        return GpuConfig(vendor, arch, described, cxx_flags, link_flags, status)

    return GpuConfig(status="none: " + ("; ".join(reasons) or "nothing tried"))


def mkdirs(*dirs):
    for directory in dirs:
        directory.mkdir(exist_ok=True, parents=True)


def main():
    mkdirs(DEF_DIR, DOT_PHARE_DIR, GPU_DIR)
    mpi_results = config_mpi()
    gpu = config_gpu()
    gen_system_file(gpu)
    write_local_cmake_file(mpi_results, gpu)


if __name__ == "__main__":
    main()
