// Unrelated module with its own TLS, dlopened after the bench module to bump the glibc DTV
//  generation ("TLS access slows after an unrelated dlopen", python dlopens many modules).

#include <cstdint>

thread_local std::uint64_t tls_dummy_value = 0;

extern "C" __attribute__((visibility("default"))) std::uint64_t tls_dummy_touch()
{
    return ++tls_dummy_value;
}
