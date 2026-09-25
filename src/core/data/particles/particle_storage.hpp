#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP


#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ostream>


// Reduced precision storage types for particle attributes.
//  Values are decoded to double on read and encoded on assignment, so all arithmetic
//  is done in double and only the stored representation loses precision.
//  Storage is a byte array (alignof == 1), so Bytes need not be a power of two.
//  Encoding rounds to nearest, and encode(decode(x)) == x bitwise.

namespace PHARE::core
{
static_assert(std::endian::native == std::endian::little,
              "particle storage types assume little endian");


template<typename Derived>
struct ParticleStorageOps
{
    template<typename T>
    Derived& operator+=(T const& v)
    {
        auto& self  = static_cast<Derived&>(*this);
        return self = static_cast<double>(self) + static_cast<double>(v);
    }

    template<typename T>
    Derived& operator-=(T const& v)
    {
        auto& self  = static_cast<Derived&>(*this);
        return self = static_cast<double>(self) - static_cast<double>(v);
    }

    template<typename T>
    Derived& operator*=(T const& v)
    {
        auto& self  = static_cast<Derived&>(*this);
        return self = static_cast<double>(self) * static_cast<double>(v);
    }

    friend std::ostream& operator<<(std::ostream& out, Derived const& d)
    {
        return out << static_cast<double>(d);
    }
};


/** The Bytes most significant bytes of an IEEE-754 double:
 *   sign, 11 exponent bits and 8 * Bytes - 12 mantissa bits.
 *   Dropped mantissa bits are rounded to nearest, ties to even.
 */
template<std::size_t Bytes>
struct TruncatedDouble : ParticleStorageOps<TruncatedDouble<Bytes>>
{
    static_assert(Bytes >= 2 and Bytes <= 7);
    static constexpr std::size_t bytes = Bytes;
    static constexpr int dropped_bits  = 64 - 8 * Bytes;

    constexpr TruncatedDouble() = default;
    TruncatedDouble(double const v) { *this = v; }

    TruncatedDouble& operator=(double const v)
    {
        auto u                  = std::bit_cast<std::uint64_t>(v);
        std::uint64_t const lsb = (u >> dropped_bits) & 1;
        // a carry out of the mantissa increments the exponent, which is the correct rounding
        u += (std::uint64_t{1} << (dropped_bits - 1)) - 1 + lsb;
        u >>= dropped_bits;
        std::memcpy(data.data(), &u, Bytes);
        return *this;
    }

    operator double() const
    {
        std::uint64_t u = 0;
        std::memcpy(&u, data.data(), Bytes);
        return std::bit_cast<double>(u << dropped_bits);
    }

    std::array<std::uint8_t, Bytes> data{};
};


/** Unsigned fixed point for values in [0, 1), resolution 2^-(8 * Bytes).
 *   Rounds to nearest and clamps to [0, 1 - 2^-(8 * Bytes)], so a value that would round
 *   up to 1 can never wrap to 0 without its cell being incremented.
 *   Bytes <= 6 keeps the decoded value exactly representable as a double.
 */
template<std::size_t Bytes>
struct FixedPointUnit : ParticleStorageOps<FixedPointUnit<Bytes>>
{
    static_assert(Bytes >= 1 and Bytes <= 6);
    static constexpr std::size_t bytes = Bytes;
    static constexpr int bits          = 8 * Bytes;
    static constexpr std::uint64_t max = (std::uint64_t{1} << bits) - 1;
    static constexpr double scale      = static_cast<double>(std::uint64_t{1} << bits);

    constexpr FixedPointUnit() = default;
    FixedPointUnit(double const v) { *this = v; }

    FixedPointUnit& operator=(double const v)
    {
        double const s        = std::nearbyint(v * scale);
        std::uint64_t const u = !(s > 0)                        ? 0 // also NaN
                                : s >= static_cast<double>(max) ? max
                                                                : static_cast<std::uint64_t>(s);
        std::memcpy(data.data(), &u, Bytes);
        return *this;
    }

    operator double() const
    {
        std::uint64_t u = 0;
        std::memcpy(&u, data.data(), Bytes);
        return static_cast<double>(u) / scale;
    }

    std::array<std::uint8_t, Bytes> data{};
};



// Bytes == 8 is plain double
template<std::size_t Bytes>
struct particle_delta_storage
{
    using type = FixedPointUnit<Bytes>;
};
template<>
struct particle_delta_storage<8>
{
    using type = double;
};

template<std::size_t Bytes>
struct particle_velocity_storage
{
    using type = TruncatedDouble<Bytes>;
};
template<>
struct particle_velocity_storage<8>
{
    using type = double;
};


#ifndef PHARE_PARTICLE_DELTA_BYTES
#define PHARE_PARTICLE_DELTA_BYTES 6
#endif

#ifndef PHARE_PARTICLE_V_BYTES
#define PHARE_PARTICLE_V_BYTES 6
#endif

using ParticleDelta_t = typename particle_delta_storage<PHARE_PARTICLE_DELTA_BYTES>::type;
using ParticleV_t     = typename particle_velocity_storage<PHARE_PARTICLE_V_BYTES>::type;


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP */
