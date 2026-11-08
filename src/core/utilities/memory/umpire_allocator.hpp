#ifndef PHARE_CORE_UTILITIES_MEMORY_UMPIRE_ALLOCATOR_HPP
#define PHARE_CORE_UTILITIES_MEMORY_UMPIRE_ALLOCATOR_HPP

#include "core/def/phare_config.hpp"

#if PHARE_HAVE_UMPIRE

#include <new>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace PHARE
{

// stateless (default constructible) std:: allocator over umpire's "PHARE::data_allocator"
//  construct = false: default construct is a no-op (like KokkosAllocator / mkn NoConstructAllocator),
//   elements are not constructed on host by resize, for particles. explicit construction still happens
//  construct = true: std:: semantics, resize value-initializes (e.g. zeroed counters), for fields
template<typename T, bool construct_ = true>
class UmpireAllocator
{
    using This = UmpireAllocator<T, construct_>;

    static auto& allocator()
    {
        static auto alloc
            = umpire::ResourceManager::getInstance().getAllocator("PHARE::data_allocator");
        return alloc;
    }

public:
    using pointer         = T*;
    using reference       = T&;
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;

    UmpireAllocator() = default;
    template<typename U>
    UmpireAllocator(UmpireAllocator<U, construct_> const&) noexcept
    {
    }

    template<typename U>
    struct rebind
    {
        using other = UmpireAllocator<U, construct_>;
    };

    T* allocate(std::size_t const n) const
    {
        if (n == 0)
            return nullptr;

        void* ptr = allocator().allocate(sizeof(T) * n);

        if (!ptr)
            throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            allocator().deallocate(p);
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
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value)
    {
        if constexpr (construct_)
            ::new ((void*)ptr) U();
    }

    bool operator!=(This const& that) const { return !(*this == that); }

    bool operator==(This const& /*that*/) const
    {
        return true; // stateless
    }
};

} // namespace PHARE

#endif /* PHARE_HAVE_UMPIRE */

#endif /* PHARE_CORE_UTILITIES_MEMORY_UMPIRE_ALLOCATOR_HPP */
