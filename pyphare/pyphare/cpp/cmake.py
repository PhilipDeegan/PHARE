def decode_bytes(input):
    return input.decode("ascii", errors="ignore")


def run(cmd, shell=True, capture_output=True, check=True, print_cmd=True):
    """
    https://docs.python.org/3/library/subprocess.html
    """
    import subprocess

    if print_cmd:
        print(f"running: {cmd}")
    try:
        return subprocess.run(
            cmd, shell=shell, capture_output=capture_output, check=check
        )
    except subprocess.CalledProcessError as e:  # only triggers on failure if check=True
        raise RuntimeError(decode_bytes(e.stderr)) from e


def build(*args):
    import os

    # modules are several TUs (see src/simulator/simulator_decl.hpp): an unbounded -j can OOM
    jobs = os.environ.get("THREADS", "4")
    cmd = f"cmake --build build -j {jobs} -- " + " ".join(args)
    run(cmd, capture_output=False)
