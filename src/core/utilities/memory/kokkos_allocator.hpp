#ifndef PHARE_CORE_UTILITIES_MEMORY_KOKKOS_ALLOCATOR_HPP
#define PHARE_CORE_UTILITIES_MEMORY_KOKKOS_ALLOCATOR_HPP

#include "core/def/kokkos_tools.hpp"
#include "core/def/phare_config.hpp"

#if PHARE_HAVE_KOKKOS

#include <new>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace PHARE
{

template<typename T>
class KokkosAllocator
{
    using This = KokkosAllocator<T>;

public:
    using pointer         = T*;
    using reference       = T&;
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind
    {
        using other = KokkosAllocator<U>;
    };

    T* allocate(std::size_t const n) const
    {
        if (n == 0)
            return nullptr;

        void* ptr = Kokkos::kokkos_malloc(sizeof(T) * n);

        if (!ptr)
            throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            Kokkos::kokkos_free(p);
    }
    void deallocate(T* const p, std::size_t /*n*/) noexcept
    { // needed from std::
        deallocate(p);
    }

    template<typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        ::new ((void*)ptr) U(std::forward<Args>(args)...);
    }
    template<typename U>
    void construct(U* /*ptr*/) noexcept(std::is_nothrow_default_constructible<U>::value)
    {
    }

    bool operator!=(This const& that) const { return !(*this == that); }

    bool operator==(This const& /*that*/) const
    {
        return true; // stateless
    }
};

} // namespace PHARE

#endif /* PHARE_HAVE_KOKKOS */

#endif /* PHARE_CORE_UTILITIES_MEMORY_KOKKOS_ALLOCATOR_HPP */
