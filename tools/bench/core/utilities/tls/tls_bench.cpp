// Driver: dlopens each module built by run.sh (RTLD_NOW | RTLD_LOCAL, as python), times every
//  tls_run_<kind> in ns per access, on the main thread and on a fresh thread, then again after
//  dlopening an unrelated TLS module (tls_dummy) without touching its TLS.
//
//  usage: tls_bench [--calls N] [--reps R] [--dummy path/to/libtls_dummy.so] module.so...

#include <dlfcn.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

namespace
{
using run_fn = std::uint64_t (*)(std::uint64_t);

char const* const kinds[] = {"global",      "block_obj", "scratch_ptr", "scratch_ptr_ie",
                             "ns_obj",      "ns_ptr",    "static_ptr"};

struct Result
{
    double min_ns, med_ns;
};

std::uint64_t volatile sink = 0;

Result time_it(run_fn const fn, std::uint64_t const calls, std::size_t const reps)
{
    sink = sink + fn(1000); // warmup + first access (lazy TLS block allocation / init) untimed

    std::vector<double> ns(reps);
    for (auto& t : ns)
    {
        auto const t0 = std::chrono::steady_clock::now();
        sink          = sink + fn(calls);
        auto const t1 = std::chrono::steady_clock::now();
        t = std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(calls);
    }
    std::sort(ns.begin(), ns.end());
    return {ns.front(), ns[reps / 2]};
}

Result time_on_thread(run_fn const fn, std::uint64_t const calls, std::size_t const reps)
{
    Result r{};
    std::thread{[&] { r = time_it(fn, calls, reps); }}.join();
    return r;
}

std::string basename_of(std::string const& path)
{
    auto const pos = path.find_last_of('/');
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

void print(std::string const& mod, char const* kind, char const* phase, Result const r)
{
    std::printf("%-34s %-15s %-14s %8.3f %8.3f\n", mod.c_str(), kind, phase, r.min_ns, r.med_ns);
}

} // namespace

int main(int argc, char** argv)
{
    std::uint64_t calls = 100'000'000;
    std::size_t reps    = 7;
    char const* dummy   = nullptr;
    std::vector<std::string> modules;

    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--calls") && i + 1 < argc)
            calls = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--reps") && i + 1 < argc)
            reps = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--dummy") && i + 1 < argc)
            dummy = argv[++i];
        else
            modules.emplace_back(argv[i]);
    }
    if (modules.empty() || reps == 0 || calls == 0)
    {
        std::fprintf(stderr, "usage: %s [--calls N] [--reps R] [--dummy lib.so] module.so...\n",
                     argv[0]);
        return 1;
    }

    std::printf("# calls %llu reps %zu, ns per access (min, median)\n",
                static_cast<unsigned long long>(calls), reps);
    std::printf("%-34s %-15s %-14s %8s %8s\n", "module", "kind", "phase", "min", "median");

    // each module in its own phase sequence; the dummy dlopen is per module so every module's
    //  "+dlopen" rows see exactly one generation bump after its own load
    for (auto const& path : modules)
    {
        void* const handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!handle)
        {
            std::fprintf(stderr, "dlopen %s: %s\n", path.c_str(), dlerror());
            return 1;
        }
        auto const mod = basename_of(path);

        std::vector<std::pair<char const*, run_fn>> fns;
        for (auto const* kind : kinds)
            if (auto* sym = dlsym(handle, (std::string{"tls_run_"} + kind).c_str()))
                fns.emplace_back(kind, reinterpret_cast<run_fn>(sym));

        for (auto const& [kind, fn] : fns)
        {
            print(mod, kind, "main", time_it(fn, calls, reps));
            print(mod, kind, "thread", time_on_thread(fn, calls, reps));
        }

        if (dummy)
        {
            void* const dh = dlopen(dummy, RTLD_NOW | RTLD_LOCAL);
            if (!dh)
            {
                std::fprintf(stderr, "dlopen %s: %s\n", dummy, dlerror());
                return 1;
            }
            for (auto const& [kind, fn] : fns)
            {
                print(mod, kind, "main+dlopen", time_it(fn, calls, reps));
                print(mod, kind, "thread+dlopen", time_on_thread(fn, calls, reps));
            }
            dlclose(dh);
        }
        // modules stay loaded: dlclose of a module with live thread_local objects is a mess
    }
    return static_cast<int>(sink & 0); // keep sink alive
}
