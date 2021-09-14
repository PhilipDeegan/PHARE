#ifndef PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H

#ifndef HAVE_UMPIRE
#error // expected
#endif

#include <atomic>
#include <cstddef>

#include "kul/gpu.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"

namespace PHARE::core::llnl
{
template<typename Particle>
struct ABufferedParticleVector
{
    auto get_allocator()
    {
        auto& rm = umpire::ResourceManager::getInstance();
        assert(rm.isAllocator("PHARE::data_allocator"));
        return rm.getAllocator("PHARE::data_allocator");
    }


    ABufferedParticleVector(std::size_t size_, double buffer_by_, double realloc_by_)
        : size{size_}
        , buffer_by{realloc_by_}
        , realloc_by{realloc_by_}
        , capacity{static_cast<std::size_t>(size + size * buffer_by)}
        , allocator_{get_allocator()}
        , particles{allocator_}
    {
        particles.resize(capacity);
        assert(particles.size() == capacity);
        assert(buffer_by < 1 and buffer_by > 0);
        assert(realloc_by < 1 and realloc_by > 0);

        // info = std::make_unique<kul::gpu::DeviceMem<std::size_t>>(std::vector(1, size));
        // //KLOG(INF) << info->p;
        // //KLOG(INF) << info.get();
    }

    bool check() __host__
    {
        auto n_elements = (*info)()[0];

        bool realloc_more = n_elements >= size + size * realloc_by;
        if (realloc_more)
        {
            size = n_elements;
            capacity += size * buffer_by;
            particles.reserve(capacity);
        }

        // bool realloc_less = false;
        // if (realloc_less)
        // {
        //     // TODO?
        // }
    }

    std::size_t size;
    double buffer_by, realloc_by;
    std::size_t capacity;
    std::unique_ptr<kul::gpu::DeviceMem<std::size_t>> info;

    umpire::TypedAllocator<double> allocator_;
    std::vector<double, umpire::TypedAllocator<double>> particles;
};


template<typename Particle>
struct BufferedParticleVector : ABufferedParticleVector<Particle>
{
    using Super    = ABufferedParticleVector<Particle>;
    using iterator = Particle*;
    using Super::info;
    using Super::particles;

    BufferedParticleVector(std::size_t size, double buffer_by = .1,
                           double realloc_by = .05) __host__
        : ABufferedParticleVector<Particle>{size, buffer_by, realloc_by}
    {
    }

    std::size_t size() const __host__ { return Super::size; }
    std::size_t size() const __device__ { return info[0]; }

    auto& operator[](std::size_t i) const __device__ { return particles[i]; }
    auto& operator[](std::size_t i) __device__ { return particles[i]; }

    void push_back(Particle const& particle) __device__
    {
        assert(false);
        // auto index       = atomicAdd(&info[0], 1);
        // particles[index] = particle;
    }
    void push_back(Particle&& particle) __device__ { push_back(particle); }

    void erase(std::size_t index) __device__
    {
        assert(false);
        // auto top = atomicSub(&info[0], 1);
        // if (index != top)
        //     particles[index] = top;
    }
};

template<typename Particle, typename Vector_ = BufferedParticleVector<Particle>>
class ParticleArray
{
public:
    static constexpr bool is_host_mem   = false;
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = Particle::dimension;
    using Particle_t                    = Particle;
    using Vector                        = Vector_;
    using iterator                      = typename Vector::iterator;
    using value_type                    = Particle_t;

    ParticleArray() { std::cout << __FILE__ << " " << __LINE__ << std::endl; };
    ParticleArray(std::size_t size)
        : vector{std::make_unique<Vector>(size)}
    {
        std::cout << __FILE__ << " " << __LINE__ << " " << size << std::endl;
    }

    std::size_t size() const
    {
        check();
        return vector->size();
    }

    void clear() { vector.release(); }

    auto& operator=(std::vector<Particle_t>&& input)
    {
        // this->particles = std::move(vector);
        vector = std::make_unique<Vector>(input.size());
        std::cout << __FILE__ << " " << __LINE__ << " " << input.size() << std::endl;

#if defined(HAVE_RAJA)
        RAJA::resources::Cuda{}.memcpy(
            /*device pointer*/ vector->particles.data(),
            /*host pointer*/ input.data(),
            /*size in bytes*/ sizeof(Particle_t) * input.size());
#else
        static_assert(false, "no other impl defined");
#endif
        // KLOG(INF);
        input.clear(); // probably dies here but just in case;
        return *this;
    }

    operator bool() const { return vector != nullptr; }
    void check() const { assert(bool{*this}); }

    // auto& operator[](std::size_t i) _PHARE_ALL_FN_
    // {
    //     check();
    //     return vector[i];
    // }
    // auto& operator[](std::size_t i) const _PHARE_ALL_FN_
    // {
    //     check();
    //     return vector[i];
    // }

    auto data() const { return vector->particles.data(); }
    auto data() { return vector->particles.data(); }
    void erase(std::size_t index)
    {
        assert(false);
        // vector->erase(index);
    }

    void push_back(Particle&& p) { vector.push_back(p); }
    void push_back(Particle const& p) { vector.push_back(p); }

    void swap(ParticleArray<Particle>& that)
    {
        this->vector->particles.swap(that.vector->particles);
    }

    auto operator()()()
    {
        check();
        std::vector<Particle> particles(size());

#if defined(HAVE_RAJA)
        RAJA::resources::Cuda{}.memcpy( // expects CUDA unified mem
            /*dst host pointer*/ particles.data(),
            /*src dev pointer*/ vector->particles.data(),
            /*size in bytes*/ sizeof(Particle) * size());
#else
        static_assert(false, "no other impl defined");
#endif

        return std::move(particles);
    }

private:
    std::unique_ptr<Vector> vector;
};

} // namespace PHARE::core::llnl


namespace PHARE::core
{
template<typename Particle>
void empty(llnl::ParticleArray<Particle>& array)
{
    assert(false);
    array.clear();
}

template<typename Particle>
void swap(llnl::ParticleArray<Particle>& array1, llnl::ParticleArray<Particle>& array2)
{
    assert(false);
    array1.swap(array2);
}

template<typename Particle>
void append(llnl::ParticleArray<Particle> const& src, llnl::ParticleArray<Particle>& dst)
{
    assert(false);
}

} // namespace PHARE::core

#endif /*PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H*/
