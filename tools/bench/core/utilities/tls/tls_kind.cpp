// One access pattern per object file (-DTLS_KIND=<n>), all linked into one -fPIC -shared module
//  that the driver dlopens, like python dlopens our cpp_* modules.
//  Separate TUs on purpose: with clang one dynamically initialised thread_local in a TU makes
//  every thread_local of that TU go through init calls, which would contaminate the pointer kinds.
//
//  Each kind exports `extern "C" std::uint64_t tls_run_<name>(std::uint64_t n)`: n calls of a
//  noinline accessor, i.e. one TLS lookup per call (the "looked up once per box/tile" case),
//  not a hoisted loop.

#include <cstdint>
#include <vector>

#ifndef TLS_KIND
#error "define TLS_KIND"
#endif

#define TLS_CAT_(a, b) a##b
#define TLS_CAT(a, b) TLS_CAT_(a, b)
#define TLS_NS TLS_CAT(kind_, TLS_KIND) // distinct per kind: no ODR merging across objects

namespace tls_bench::TLS_NS
{
struct Scratch // non trivial, like the refiner / ion_updater scratch vectors
{
    std::vector<double> buf;
    std::uint64_t n = 0;
};

#define TLS_NOINLINE __attribute__((noinline))

#if TLS_KIND == 0 // baseline: plain global, no TLS (call + increment cost)
#define TLS_NAME global
Scratch g_scratch;
TLS_NOINLINE Scratch& access()
{
    return g_scratch;
}

#elif TLS_KIND == 1 // (a) block-scope thread_local std::vector-holding object (guard per access)
#define TLS_NAME block_obj
TLS_NOINLINE Scratch& access()
{
    thread_local Scratch s;
    return s;
}

#elif TLS_KIND == 2 // (b) proposed thread_scratch<T>(): constant-init thread_local T*, lazy new
#define TLS_NAME scratch_ptr
template<typename T>
T& thread_scratch()
{
    thread_local T* p = nullptr;
    if (!p)
        p = new T{}; // once per thread, never freed
    return *p;
}
TLS_NOINLINE Scratch& access()
{
    return thread_scratch<Scratch>();
}

#elif TLS_KIND == 3 // (c) (b) + per-variable initial-exec
#define TLS_NAME scratch_ptr_ie
template<typename T>
T& thread_scratch()
{
    [[gnu::tls_model("initial-exec")]] thread_local T* p = nullptr;
    if (!p)
        p = new T{};
    return *p;
}
TLS_NOINLINE Scratch& access()
{
    return thread_scratch<Scratch>();
}

#elif TLS_KIND == 4 // namespace-scope external thread_local object (TLS wrapper/init path)
#define TLS_NAME ns_obj
thread_local Scratch tl_ns_obj;
TLS_NOINLINE Scratch& access()
{
    return tl_ns_obj;
}

#elif TLS_KIND == 5 // namespace-scope external constant-init thread_local pointer
#define TLS_NAME ns_ptr
thread_local Scratch* tl_ns_ptr = nullptr;
TLS_NOINLINE Scratch& access()
{
    if (!tl_ns_ptr)
        tl_ns_ptr = new Scratch{};
    return *tl_ns_ptr;
}

#elif TLS_KIND == 6 // file-local (internal linkage) constant-init thread_local pointer
#define TLS_NAME static_ptr
namespace
{
    thread_local Scratch* tl_static_ptr = nullptr;
}
TLS_NOINLINE Scratch& access()
{
    if (!tl_static_ptr)
        tl_static_ptr = new Scratch{};
    return *tl_static_ptr;
}

#else
#error "unknown TLS_KIND"
#endif

} // namespace tls_bench::TLS_NS

extern "C" __attribute__((visibility("default"))) std::uint64_t
TLS_CAT(tls_run_, TLS_NAME)(std::uint64_t const n)
{
    std::uint64_t sum = 0;
    for (std::uint64_t i = 0; i < n; ++i)
    {
        auto& s = tls_bench::TLS_NS::access();
        sum += ++s.n;
        asm volatile("" ::: "memory"); // no merging / hoisting of calls across iterations
    }
    return sum;
}
