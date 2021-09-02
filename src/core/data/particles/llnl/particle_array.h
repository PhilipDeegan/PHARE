#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H

#include <atomic>
#include <cstddef>

#include "kul/gpu.hpp"

namespace PHARE::core::llnl
{
template<bool GPU, typename Particle>
struct ABufferedParticleVector
{
    ABufferedParticleVector() = delete;
};

template<typename Particle>
struct ABufferedParticleVector<false, Particle>
{
    ABufferedParticleVector(std::size_t size_, double buffer_by_, double realloc_by_)
        : size{size_}
        , buffer_by{realloc_by_}
        , realloc_by{realloc_by_}
        , capacity{size + size * buffer_by}
    {
        assert(buffer_by < 1 and buffer_by > 0);
        assert(realloc_by < 1 and realloc_by > 0);

        info      = std::make_unique<kul::gpu::DeviceMem<std::size_t>>(1, size);
        particles = std::make_unique<kul::gpu::DeviceMem<Particle>>(v_info[3]);
    }

    bool check() __host__
    {
        auto n_elements = info()[0];

        bool realloc_more = n_elements >= size + size * realloc_by;
        if (realloc_more)
        {
            size = n_elements;
            capacity += size * buffer_by;
            auto new_buffer = std::make_unique<kul::gpu::DeviceMem<Particle>>(capacity);
            // TODO copy old to new in kernel?
            particles = std::move(new_buffer);
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
    std::unique_ptr<kul::gpu::DeviceMem<Particle>> particles;
};
template<typename Particle>
struct ABufferedParticleVector<true, Particle>
{
    std::size_t* info;
    Particle* particles;
};

template<typename Particle, bool GPU = false>
struct BufferedParticleVector : ABufferedParticleVector<GPU, Particle>
{
    using Super    = ABufferedParticleVector<GPU, Particle>;
    using iterator = Particle*;
    using gpu_t    = BufferedParticleVector<Particle, true>;
    using Super::info;
    using Super::particles;

    BufferedParticleVector(std::size_t size, double buffer_by = .1,
                           double realloc_by = .05) __host__
        : ABufferedParticleVector<false, Particle>{size, buffer_by}
    {
    }

    auto operator()() __host__ { return Super::template alloc<gpu_t>(*info, *particles); }

    std::size_t size() const __host__ { return Super::size; }
    std::size_t size() const __device__ { return info[0]; }

    auto& operator[](std::size_t i) const __device__ { return particles[i]; }
    auto& operator[](std::size_t i) __device__ { return particles[i]; }


    void push_back(Particle const& particle) __device__
    {
        auto index       = atomicAdd(&info[0], 1);
        particles[index] = particle;
    }

    void push_back(Particle&& particle) __device__
    {
        auto index       = atomicAdd(&info[0], 1);
        particles[index] = particle;
    }

    void erase(std::size_t index) __device__
    {
        auto top = atomicSub(&info[0], 1);
        if (index != top)
            particles[index] = top;
    }
};

template<typename Particle, typename Vector_ = BufferedParticleVector>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = Particle::dimension;
    using Particle_t                    = Particle;
    using Vector                        = Vector_;
    using iterator                      = typename Vector::iterator;
    using value_type                    = Particle_t;

    ParticleArray() = delete;
    ParticleArray(std::size_t size)
        : particles(size)
    {
    }

    std::size_t size() const { return vector.size(); }

    auto& operator[](std::size_t i) { return vector[i]; }
    auto& operator[](std::size_t i) const { return vector[i]; }

    void erase(std::size_t index) { return vector.erase(index); }

    void push_back(Particle_t&& p) { vector.push_back(p); }
    void push_back(Particle_t const& p) { vector.push_back(p); }

    void swap(ParticleArray<dim>& that) { this->vector.particles.swap(that.vector.particles); }

private:
    Vector vector;
};

} // namespace PHARE::core::llnl


namespace PHARE::core
{
template<std::size_t dim>
void empty(llnl::ParticleArray<dim>& array)
{
    static_assert(false);
    array.clear();
}

template<std::size_t dim>
void swap(llnl::ParticleArray<dim>& array1, llnl::ParticleArray<dim>& array2)
{
    static_assert(false);
    array1.swap(array2);
}

} // namespace PHARE::core

#endif

