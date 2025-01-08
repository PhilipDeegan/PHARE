#ifndef PHARE_CORE_VECTOR_HPP
#define PHARE_CORE_VECTOR_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/memory.hpp"

#include <vector>

namespace PHARE
{


auto constexpr default_allocator_mode()
{
    if constexpr (CompileOptions::WithUmpire or CompileOptions::WithMknGpu)
        return AllocatorMode::GPU_UNIFIED;
    return AllocatorMode::CPU;
}

template<auto allocator_mode>
bool constexpr allocator_mode_supported()
{
    if constexpr (allocator_mode == AllocatorMode::CPU)
        return true;
    else if constexpr (/*allocator_mode == AllocatorMode::GPU
                       or */
                       allocator_mode == AllocatorMode::GPU_UNIFIED)
        return CompileOptions::WithMknGpu or CompileOptions::WithUmpire;
    return false;
}

template<auto _allocator_mode, typename Actual>
struct Allocator
{
    auto static constexpr mode = _allocator_mode;
    using value_type           = Actual;

    auto& operator()() { return actual; }

    Actual actual{};
};

template<typename Type, auto allocator_mode, std::uint8_t mode = 0>
auto constexpr allocator()
{
    static_assert(std::is_same_v<decltype(allocator_mode), AllocatorMode>);
    static_assert(allocator_mode_supported<allocator_mode>());

    if constexpr (allocator_mode == AllocatorMode::CPU)
        return Allocator<allocator_mode, typename std::vector<Type>::allocator_type>{};
    else if constexpr (/*allocator_mode == AllocatorMode::GPU
                       or */
                       allocator_mode == AllocatorMode::GPU_UNIFIED)
    {
        if constexpr (CompileOptions::WithMknGpu)
        {
            PHARE_WITH_MKN_GPU({
                if constexpr (mode == 0) // for fields
                    return Allocator<allocator_mode, mkn::gpu::ManagedAllocator<Type>>{};
                else // for particles
                    return Allocator<allocator_mode, mkn::gpu::NoConstructAllocator<Type>>{};
            });
        }
        if constexpr (CompileOptions::WithUmpire)
        {
            PHARE_WITH_UMPIRE(return Allocator<allocator_mode, umpire::TypedAllocator<Type>>{
                umpire::ResourceManager::getInstance().getAllocator("PHARE::data_allocator")});
        }

        // or compile error
    }
    throw std::runtime_error("NOOO");
}

template<typename Type, auto allocator_mode, std::uint8_t mode>
using allocator_type_t = decltype(allocator<Type, allocator_mode, mode>());

// ^ bi directionality v

template<typename Type, typename allocator_type, std::uint8_t mode = 0>
auto constexpr get_allocator_mode()
{
    if constexpr (std::is_same_v<allocator_type, typename allocator_type_t<Type, AllocatorMode::CPU,
                                                                           mode>::value_type>)
        return AllocatorMode::CPU;
    else
    {
        return AllocatorMode::GPU_UNIFIED; // do better
    }
}


template<typename Type, auto allocator_mode_ = AllocatorMode::CPU, std::uint8_t mode = 0,
         typename vector_t_
         = std::vector<Type, typename allocator_type_t<Type, allocator_mode_, mode>::value_type>>
struct Vector
{
    auto static constexpr allocator_mode = allocator_mode_;

    using value_type     = Type;
    using Allocator_t    = allocator_type_t<Type, allocator_mode, mode>;
    using allocator_type = typename Allocator_t::value_type;
    using vector_t       = vector_t_;
    static_assert(allocator_mode == get_allocator_mode<Type, allocator_type, mode>());

    Vector() = delete; // not meant to be instantiated.

    auto static make(std::size_t size, bool reserve = false)
    {
        static_assert(core::is_std_vector_v<vector_t>);
        vector_t vec; //{Allocator_t{}()};
        if (size)
        {
            if (reserve)
                vec.reserve(size);
            else
                vec.resize(size);
        }
        return vec;
    }


    auto static make(std::size_t size, Type const& val)
    {
        auto vec = make(size);
        fill(vec, val);
        return vec;
    }


    template<typename Alloc0>
    void static copy(vector_t& dst, std::vector<Type, Alloc0> const& src)
    {
        auto static constexpr src_mode = get_allocator_mode<Type, Alloc0>();

        if (dst.size() != src.size())
            dst.resize(src.size());

        if constexpr (allocator_mode == src_mode)
            dst = src;
        else
            PHARE::core::gpu::copy(dst.data(), src.data(), src.size());
    }


    template<typename Alloc0>
    void static move(vector_t& dst, std::vector<Type, Alloc0>& src)
    {
        auto static constexpr src_mode = get_allocator_mode<Type, Alloc0>();

        if constexpr (allocator_mode == AllocatorMode::CPU and src_mode == AllocatorMode::CPU)
            dst = std::move(src);

        else
            copy(dst, src);
    }

    template<typename Type0, typename Alloc0>
    auto static from(std::vector<Type0, Alloc0> const& src)
    {
        static_assert(std::is_same_v<Type, Type0>);

        auto static constexpr src_mode = get_allocator_mode<Type0, Alloc0>();


        if constexpr (src_mode == AllocatorMode::CPU)
        {
            return src;
        }
        else
        {
            // allocations shouldn't happen on GPU, so assume CPU
            auto dst = make(src.size());
            copy(dst, src);
            return dst;
        }
    }

    template<typename Type0, typename Alloc0>
    auto static from(std::vector<Type0, Alloc0>&& src)
    {
        static_assert(std::is_same_v<Type, Type0>);
        auto static constexpr src_mode = get_allocator_mode<Type0, Alloc0>();

        if constexpr (allocator_mode == AllocatorMode::CPU and src_mode == AllocatorMode::CPU)
            return std::move(src);

        else
            return from(src);
    }


    auto static fill(vector_t& vec, Type const& val)
    {
        if constexpr (allocator_mode == AllocatorMode::CPU)
            std::fill(vec.begin(), vec.end(), val);

        else
            PHARE::core::gpu::fill(vec, val);
    }
};


template<typename T, auto alloc_mode_ = AllocatorMode::CPU, std::uint8_t mode = 0>
struct MinimizingVector
{
    using Allocator_t    = allocator_type_t<T, alloc_mode_, mode>;
    using allocator_type = typename Allocator_t::value_type;

    template<bool copy_old = true>
    auto& get(std::size_t const& s)
    {
        if (s < v.capacity() * percentile)
            ++_c;
        else
            _c = 0;

        if (_c == period)
        {
            std::vector<T> r(v.capacity() * realloc_to);
            if constexpr (copy_old)
                r = v;
            v  = std::move(r);
            _c = 0;
        }

        v.resize(s);
        return v;
    }

    double const percentile  = .50;
    double const realloc_to  = .60;
    std::size_t const period = 100;

    allocator_type alloc{};
    std::vector<T> v{alloc};
    std::uint16_t _c = 0;
};


} // namespace PHARE


#endif /* PHARE_CORE_VECTOR_HPP */
