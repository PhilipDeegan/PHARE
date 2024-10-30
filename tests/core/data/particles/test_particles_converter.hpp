#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP

#include "phare_core.hpp"

#include <string>
#include <cassert>
#include <functional>

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "core/def/detail/mkn_gpu.hpp"

namespace PHARE::core
{

template<auto type, typename Src, typename Dst>
struct ParticlesConverter
{
    static_assert(Src::dimension == Dst::dimension);
    auto constexpr static src_alloc_mode = Src::alloc_mode;
    auto constexpr static dst_alloc_mode = Dst::alloc_mode;

    auto constexpr static src_layout_mode = Src::layout_mode;
    auto constexpr static dst_layout_mode = Dst::layout_mode;


    template<typename GridLayout>
    auto operator()(Src const& src, GridLayout const& layout)
    {
        //
        auto const ppc = [&]() { // assume constant for the moment
            for (auto const& bix : src.box())
                if (auto const n = src.nbr_particles_in(*bix); n > 0)
                    return n;
            throw std::runtime_error("ppc 0");
        }();
        PHARE_LOG_LINE_SS(ppc);

        auto dst = make_particles<Dst>(layout);

        if constexpr (src_layout_mode == LayoutMode::AoSMapped
                      and dst_layout_mode == LayoutMode::SoAPC)
        {
            static_assert(src_alloc_mode == AllocatorMode::CPU);
            static_assert(dst_alloc_mode == AllocatorMode::GPU_UNIFIED);

            using Tmp = ParticleArray<
                Dst::dimension,
                ParticleArrayInternals<Dst::dimension, LayoutMode::SoA, StorageMode::VECTOR,
                                       src_alloc_mode, Dst::impl>>;

            auto tmp = make_particles<Tmp>(layout);
            tmp.resize(src.size());
            ParticlePacker<Src>{src}.pack(tmp);

            assert(array_equals(src.delta(0), tmp.delta(0)));

            dst.template reserve_ppc<type>(ppc);

            std::size_t tmp_start = 0;

            auto const tmp_tuple = tmp.as_tuple();
            auto finish          = [&](auto& lix) {
                auto const size = ppc;
                auto dst_tuple  = dst(lix).as_tuple();

                for_N<std::tuple_size_v<decltype(dst_tuple)>>([&](auto vi) {
                    auto const& tmp_vec = std::get<vi>(tmp_tuple);
                    auto& vec           = std::get<vi>(dst_tuple);
                    vec.resize(size);
                    PHARE_WITH_MKN_GPU(                 //
                        mkn::gpu::copy(                 //
                            vec.data(),                 //
                            tmp_vec.data() + tmp_start, //
                            size                        //
                            )                           //
                    );

                    // PHARE_LOG_LINE_SS(vec.size());
                });
                tmp_start += size;
            };

            auto const amrbox = dst.safe_box();
            auto const lclbox = dst.local_box();
            auto amr_it       = amrbox.begin();
            auto lcl_it       = lclbox.begin();

            for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
            {
                if (!array_equals(**amr_it, tmp.iCell(tmp_start)))
                    continue;

                // PHARE_LOG_LINE_SS(*amr_it << " " << Point{*tmp_cell});

                finish(*lcl_it);
            }

            dst.template sync<2, type>();
            // PHARE_LOG_LINE_SS(dst.size() << " " << src.size() << " " << tmp.size());
            assert(dst.size() == src.size());

            bool valid = false;
            amr_it     = amrbox.begin();
            lcl_it     = lclbox.begin();

            for (; amr_it != amrbox.end(); ++amr_it, ++lcl_it)
            {
                if (!array_equals(**amr_it, tmp.iCell(0)))
                    continue;

                PHARE_LOG_LINE_SS(Point{src.delta(0)} << " " << Point{dst(*lcl_it).delta(0)});

                valid = array_equals(src.delta(0), dst(*lcl_it).delta(0));
            }

            auto const& last = [&]() -> auto& {
                static_assert(std::is_same_v<decltype(type), ParticleType>);

                if constexpr (type == ParticleType::Domain)
                    return dst(dst.local_cell(dst.box().upper));
                else if (type == ParticleType::Ghost)
                    return dst(dst.local_cell(dst.ghost_box().upper));
            }();

            assert(last.size() == ppc);
            PHARE_LOG_LINE_SS(Point{src.delta(src.size() - 1)}
                              << " " << Point{last.delta(last.size() - 1)});

            assert(valid);
            return dst;
        }


        throw std::runtime_error("unconvertible");

        return make_particles<Dst>(layout);
        // else{}
    }
};


template<auto type, typename Dst, typename Src, typename GridLayout>
auto static convert(Src const& src, GridLayout const& layout)
{
    return ParticlesConverter<type, Src, Dst>{}(src, layout);
}


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_CONVERTER_HPP */
