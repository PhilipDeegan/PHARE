#ifndef PHARE_PARTICLE_INITIALIZER_HPP
#define PHARE_PARTICLE_INITIALIZER_HPP


namespace PHARE
{
namespace core
{
    template<typename ParticleArray, typename GridLayout>
    class ParticleInitializer
    {
    public:
        virtual void loadParticles(ParticleArray& particles, GridLayout const& layout) const = 0;
        virtual ~ParticleInitializer() = default;
    };

} // namespace core
} // namespace PHARE


namespace PHARE::core
{

// if you want to handle particle init manually
template<typename ParticleArray, typename GridLayout>
class NoopParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
{
public:
    NoopParticleInitializer() {}

    void loadParticles(ParticleArray& /*particles*/, GridLayout const& /*layout*/) const override {}
};


} // namespace PHARE::core

#endif
