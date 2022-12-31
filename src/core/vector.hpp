#ifndef PHARE_CORE_VECTOR_HPP
#define PHARE_CORE_VECTOR_HPP

#include <vector>

#include "core/def.hpp"
#include "core/def/detail/umpire.hpp"
#include "core/def/detail/raja.hpp" // checks for umpire to know if gpu
#include "core/def/detail/mkn_gpu.hpp"

namespace PHARE
{
struct CompileOptions
{
    static constexpr bool WithUmpire = PHARE_HAVE_UMPIRE;
    static constexpr bool WithMknGpu = PHARE_HAVE_MKN_GPU;
    static constexpr bool WithRAJA   = PHARE_HAVE_RAJA;
};
} // namespace PHARE


namespace PHARE
{
template<typename Type>
struct Vector
{
    using value_type = Type;

#if PHARE_HAVE_MKN_GPU
    using allocator_type = mkn::gpu::ManagedAllocator<Type>;
#elif PHARE_HAVE_UMPIRE
    using allocator_type = umpire::TypedAllocator<Type>;
#else
    using allocator_type = typename std::vector<Type>::allocator_type;
#endif
    using vector_type = std::vector<Type, allocator_type>;

    template<typename Allocator>
    static auto make_allocator()
    {
        if constexpr (is_host_mem<Allocator>())
        {
            return Allocator{};
        }
        else
        {
            if constexpr (CompileOptions::WithMknGpu)
            {
                PHARE_WITH_MKN_GPU(return Allocator{});
            }
            else
            {
                PHARE_WITH_UMPIRE(return Allocator{
                    umpire::ResourceManager::getInstance().getAllocator("PHARE::data_allocator")});
            }
        }
    }

    template<typename Allocator = allocator_type>
    static auto make(std::size_t size, bool reserve = false)
    {
        vector_type vec(make_allocator<Allocator>());
        if (size)
        {
            if (reserve)
                vec.reserve(size);
            else
                vec.resize(size);
        }
        return vec;
    }


    template<typename Allocator>
    static constexpr bool is_host_mem()
    {
        return std::is_same_v<Allocator, typename std::vector<Type>::allocator_type>;
    }


    template<typename Alloc0, typename Alloc1>
    static void copy(std::vector<Type, Alloc0>& dst, std::vector<Type, Alloc1> const& src)
    {
        if (dst.size() != src.size())
            dst.resize(src.size());

        if constexpr (is_host_mem<Alloc0>() and is_host_mem<Alloc1>())
        {
            dst = src;
        }
        else
        {
            if constexpr (CompileOptions::WithRAJA)
            {
                PHARE_WITH_RAJA(PHARE::core::raja::copy(dst.data(), src.data(), src.size()));
            }
            else if (CompileOptions::WithMknGpu)
            {
                PHARE_WITH_MKN_GPU(mkn::gpu::copy(dst.data(), src.data(), src.size()));
            }
            else
                throw std::runtime_error("Vector::copy NO ALTERNATIVE");
        }
    }




    template<typename Alloc0, typename Alloc1>
    static void move(std::vector<Type, Alloc0>& dst, std::vector<Type, Alloc1>& src)
    {
        if constexpr (is_host_mem<Alloc0>() and is_host_mem<Alloc1>())
        {
            dst = std::move(src);
        }
        else
        {
            copy(dst, src);
        }
    }

    template<typename Allocator>
    static auto from(std::vector<Type, Allocator> const& that)
    {
        if constexpr (is_host_mem<Allocator>())
        {
            return that;
        }
        else
        {
            auto dst = make<Allocator>(that.size());
            copy(dst, that);
            return dst;
        }
    }

    template<typename Allocator>
    static auto from(std::vector<Type, Allocator>&& that)
    {
        if constexpr (is_host_mem<Allocator>())
        {
            auto dst = std::move(that);
            return dst;
        }
        else
        {
            return from(that);
        }
    }


    template<typename Allocator>
    static auto fill(std::vector<Type, Allocator>& vec, Type val)
    {
        if constexpr (is_host_mem<Allocator>())
        {
            std::fill(vec.begin(), vec.end(), val);
        }
        else
        {
            if constexpr (CompileOptions::WithRAJA)
            {
                PHARE_WITH_RAJA(PHARE::core::raja::set(vec.data(), val, vec.size()));
            }
            else if (CompileOptions::WithMknGpu)
            {
                PHARE_WITH_MKN_GPU(mkn::gpu::fill(vec, val));
            }
            else
                throw std::runtime_error("Vector::fill NO ALTERNATIVE");
        }
    }
};
} // namespace PHARE

#endif /* PHARE_CORE_VECTOR_HPP */
