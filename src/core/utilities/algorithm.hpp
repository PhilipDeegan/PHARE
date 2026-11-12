#ifndef PHARE_CORE_UTILITIES_ALGORITHM_HPP
#define PHARE_CORE_UTILITIES_ALGORITHM_HPP


#include "core/def.hpp"

#include "core/utilities/span.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/data/field/recentered_field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/grid/grid_tiles.hpp"

#include <array>
#include <tuple>
#include <string>
#include <algorithm>
#include <stdexcept>


namespace PHARE
{
namespace core
{
    template<std::uint32_t lhs, std::uint32_t rhs>
    NO_DISCARD constexpr std::uint32_t max()
    {
        if constexpr (lhs < rhs)
        {
            return rhs;
        }
        else if constexpr (lhs >= rhs)
        {
            return lhs;
        }
    }


    template<typename T>
    NO_DISCARD std::string to_str(T&& t)
    {
        return t.to_str();
    }



    template<typename Container, typename ContainedT = Container::value_type>
    NO_DISCARD bool notIn(ContainedT& obj, Container& list)
    {
        auto sameItem = std::find_if(std::begin(list), std::end(list), [&obj](auto& currentItem) {
            return obj->name() == currentItem->name();
        });

        return sameItem == std::end(list);
    }



} // namespace core
} // namespace PHARE


namespace PHARE::core
{

template<Spannable Span>
void average(Span const& f1, Span const& f2, Span& avg)
{
    auto const size = f1.size();
    auto const d1   = f1.data();
    auto const d2   = f2.data();
    auto av         = avg.data();
    for (std::size_t i = 0; i < size; ++i)
        av[i] = (d1[i] + d2[i]) * .5;
}

// TODO benchmark several ways of converting to primal
// https://github.com/PHAREHUB/PHARE/pull/1486#issuecomment-6401620199
template<std::size_t dim, typename PQ, typename D, auto am>
auto& convert_to_fortran_primal(      // DOES NOT WORK ON GHOST BOX!
    Field<dim, PQ, D, am>& dst,       //
    Field<dim, PQ, D, am> const& src, //
    auto const& layout                //
)
{
    bool static constexpr c_ordering = false;

    if (not all(layout.centering(dst), [](auto const c) { return c == QtyCentering::primal; }))
        throw std::runtime_error("Invalid operation, all directions must be primal");

    auto const all_primal
        = all(layout.centering(src), [](auto const c) { return c == QtyCentering::primal; });

    auto const lcl_src_box = layout.AMRToLocal(layout.AMRBoxFor(dst));
    auto dst_view          = make_array_view<c_ordering>(dst.data(), *lcl_src_box.shape());

    auto const lcl_dst_box = [&]() {
        auto box = lcl_src_box;
        box.lower -= lcl_src_box.lower; // 0
        box.upper -= lcl_src_box.lower; // shift
        return box;
    }();

    if (all_primal)
        for_each_point_zip(lcl_src_box, lcl_dst_box,
                           [&](auto const& lix, auto const& lix0) { dst_view(lix0) = src(lix); });

    else
    {
        auto const isDual = std::apply(
            [](auto... dual) { return Constexprifier{dual...}; },
            generate_from([](auto c) { return c == QtyCentering::dual; }, layout.centering(src)));

        isDual([&]<bool... dual>() {
            constexpr std::array source{(dual ? QtyCentering::dual : QtyCentering::primal)...};
            constexpr auto destination = ConstArray<QtyCentering, dim>(QtyCentering::primal);
            CompileTimeRecenteredField<source, destination, Field<dim, PQ, D, am>> const primal_src{
                src, layout};
            for_each_point_zip(lcl_src_box, lcl_dst_box, [&](auto const& lix, auto const& lix0) {
                dst_view(lix0) = primal_src(*lix);
            });
        });
    }

    return dst;
}

template<typename Field_t, typename PQ, std::size_t rank>
auto& convert_to_fortran_primal(               //
    TensorField<Field_t, PQ, rank>& dst,       //
    TensorField<Field_t, PQ, rank> const& src, //
    auto const& layout                         //
)
{
    for (std::size_t ci = 0; ci < src.size(); ++ci)
        convert_to_fortran_primal(dst[ci], src[ci], layout);
    return dst;
}

template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in, FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    std::size_t const size = in().size();
    for (std::size_t i = 0; i < size; ++i)
        std::transform(in[i]().begin(), in[i]().end(), out[i]().begin(), fn);
}


template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in0, FieldTileSet<GL, ND, PQ> const& in1,
               FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    assert(in0.isUsable());
    assert(in1.isUsable());
    assert(out.isUsable());
    auto const size = in0().size();

    for (std::size_t i = 0; i < size; ++i)
    {
        auto& i0 = in0[i]();
        auto& i1 = in1[i]();
        auto& o0 = out[i]();
        std::transform(in0[i]().begin(), in0[i]().end(), in1[i]().begin(), out[i]().begin(), fn);
    }
}


// average as jobs: (x0, x1, avg, size) over contiguous buffers, an untiled field is one buffer, a
//  tile set one per tile. Collected per field, patch or level, then run once: a device launch for
//  GPU_UNIFIED (a block per job), else on the host
namespace detail
{
    template<typename T>
    struct AverageJob
    {
        T const* x0;
        T const* x1;
        T* avg;
        std::size_t size;
    };

    // buffers are split so a whole patch isn't one block
    std::size_t constexpr static average_job_chunk = 1 << 12;

    void append_average_jobs(auto const* x0, auto const* x1, auto* avg, std::size_t const size,
                             auto& jobs)
    {
        for (std::size_t off = 0; off < size; off += average_job_chunk)
            jobs.push_back({x0 + off, x1 + off, avg + off,
                            std::min(average_job_chunk, size - off)});
    }
} // namespace detail

template<auto opts>
void average_jobs(basic::Field<opts> const& f1, basic::Field<opts> const& f2,
                  basic::Field<opts>& avg, auto& jobs)
{
    detail::append_average_jobs(f1.data(), f2.data(), avg.data(), avg.size(), jobs);
}

template<typename GL, typename ND, typename PQ>
void average_jobs(FieldTileSet<GL, ND, PQ> const& f1, FieldTileSet<GL, ND, PQ> const& f2,
                  FieldTileSet<GL, ND, PQ>& avg, auto& jobs)
{
    for (std::size_t i = 0; i < avg().size(); ++i)
        detail::append_average_jobs(f1[i]().data(), f2[i]().data(), avg[i]().data(),
                                    avg[i]().size(), jobs);
}

template<typename Field_t, std::size_t rank>
void average_jobs(basic::TensorField<Field_t, rank> const& vf1,
                  basic::TensorField<Field_t, rank> const& vf2,
                  basic::TensorField<Field_t, rank>& Vavg, auto& jobs)
{
    for (std::size_t i = 0; i < detail::tensor_field_dim_from_rank<rank>(); ++i)
        average_jobs(vf1[i], vf2[i], Vavg[i], jobs);
}

// the scalar type of a field's buffers, tiled or not
template<typename F>
auto constexpr field_scalar_type_()
{
    if constexpr (is_field_tile_set_v<F>)
        return std::type_identity<typename F::value_type::type>{};
    else
        return std::type_identity<typename F::value_type>{};
}
template<typename F>
using field_scalar_t = decltype(field_scalar_type_<std::remove_cvref_t<F>>())::type;

// a reused job list per thread, so a call doesn't allocate (managed memory for GPU_UNIFIED)
template<typename T, auto alloc_mode>
auto& average_job_list()
{
    using Job_t = detail::AverageJob<T>;
    if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
#if PHARE_HAVE_GPU
        static thread_local gpu::Vec_t<Job_t> jobs;
#else // run_average_jobs throws
        static thread_local std::vector<Job_t> jobs;
#endif
        jobs.clear();
        return jobs;
    }
    else
    {
        static thread_local std::vector<Job_t> jobs;
        jobs.clear();
        return jobs;
    }
}

template<auto alloc_mode>
void run_average_jobs(auto const& jobs)
{
    using Job_t = std::decay_t<decltype(jobs[0])>;
    if constexpr (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
#if PHARE_HAVE_GPU
        if (jobs.empty())
            return;
        using Launcher = gpu::FieldTileLauncher</*sync=*/true>;
        auto* const jp = jobs.data();
        Launcher{jobs.size()}([=] _PHARE_ALL_FN_() mutable {
            auto const& job = jp[Launcher::tile_idx()];
            Launcher::for_each_index(job.size, [&](auto const i) {
                job.avg[i] = (job.x0[i] + job.x1[i]) * .5;
            });
        });
#else
        throw std::runtime_error("average: GPU_UNIFIED needs a GPU build");
#endif
    }
    else
        for (Job_t const& job : jobs)
            for (std::size_t i = 0; i < job.size; ++i)
                job.avg[i] = (job.x0[i] + job.x1[i]) * .5;
}

template<typename GL, typename ND, typename PQ>
void average(FieldTileSet<GL, ND, PQ> const& f1, FieldTileSet<GL, ND, PQ> const& f2,
             FieldTileSet<GL, ND, PQ>& avg)
{
    using Tiles_t = FieldTileSet<GL, ND, PQ>;
    if constexpr (Tiles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        using T    = std::remove_cvref_t<decltype(*avg[0]().data())>;
        auto& jobs = average_job_list<T, Tiles_t::alloc_mode>();
        average_jobs(f1, f2, avg, jobs);
        run_average_jobs<Tiles_t::alloc_mode>(jobs);
    }
    else
    {
        core::transform(f1, f2, avg, std::plus<double>());
        core::transform(avg, avg, [](double x) { return x * 0.5; });
    }
}


template<auto opts>
void transform(basic::Field<opts> const& in, basic::Field<opts>& out, auto const fn)
{
    std::transform(in.begin(), in.end(), out.begin(), fn);
}

template<auto opts>
void transform(basic::Field<opts> const& in0, basic::Field<opts> const& in1,
               basic::Field<opts>& out, auto const fn)
{
    std::transform(in0.begin(), in0.end(), in1.begin(), out.begin(), fn);
}


template<auto opts>
void average(basic::Field<opts> const& f1, basic::Field<opts> const& f2, basic::Field<opts>& avg)
{
    if constexpr (opts.alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        using T    = std::remove_cvref_t<decltype(*avg.data())>;
        auto& jobs = average_job_list<T, opts.alloc_mode>();
        average_jobs(f1, f2, avg, jobs);
        run_average_jobs<opts.alloc_mode>(jobs);
    }
    else
    {
        core::transform(f1, f2, avg, std::plus<double>());
        core::transform(avg, avg, [](double x) { return x * 0.5; });
    }
}

template<typename Field_t, std::size_t rank>
void average(basic::TensorField<Field_t, rank> const& vf1,
             basic::TensorField<Field_t, rank> const& vf2, basic::TensorField<Field_t, rank>& Vavg)
{
    auto constexpr static N = detail::tensor_field_dim_from_rank<rank>();

    for (std::size_t i = 0; i < N; ++i)
        average(vf1[i], vf2[i], Vavg[i]);
}


template<typename Field_t>
void accumulate_field(Field_t& dst, Field_t const& src, auto const coef)
    requires(not is_field_tile_set_v<Field_t>)
{
    if (dst.size() != src.size())
        throw std::runtime_error("Cannot accumulate fields with different sizes");

    for (std::size_t i = 0; i < dst.size(); ++i)
        dst.data()[i] = src.data()[i] * coef;
}


template<typename Field_t>
void accumulate_field(Field_t& dst, Field_t const& src, auto const coef)
    requires(is_field_tile_set_v<Field_t>)
{
    auto& dst_tiles = dst();
    auto& src_tiles = src();

    if (dst_tiles.size() != src_tiles.size())
        throw std::runtime_error("Cannot accumulate fields with different number of tiles");

    for (std::size_t i = 0; i < dst_tiles.size(); ++i)
        accumulate_field(dst_tiles[i](), src_tiles[i](), coef);
}



template<typename Field_t>
void accumulate(basic::TensorField<Field_t, 1>& dst, basic::TensorField<Field_t, 1> const& src,
                auto const coeff)
{
    for (std::size_t c = 0; c < dst.size(); ++c)
        accumulate_field(dst[c], src[c], coeff);
}

template<template<typename> typename Op, typename Field_t>
void operate_field(Field_t& dst, Field_t const& src, auto&&... args)
    requires(not is_field_tile_set_v<Field_t>)
{
    using Operator = Op<typename Field_t::type>;
    for (std::size_t i = 0; i < dst.size(); ++i)
        Operator{dst.data()[i]}(src.data()[i], args...);
}

template<template<typename> typename Op, typename Field_t>
void operate_field(Field_t& dst, Field_t const& src, auto&&... args)
    requires(is_field_tile_set_v<Field_t>)
{
    auto& dst_tiles = dst();
    auto& src_tiles = src();
    if (dst_tiles.size() != src_tiles.size())
        throw std::runtime_error("operate_field: tile count mismatch");
    for (std::size_t i = 0; i < dst_tiles.size(); ++i)
        operate_field<Op>(dst_tiles[i](), src_tiles[i](), args...);
}

template<template<typename> typename Op, typename Field_t>
void operate(Field_t& dst, Field_t const& src, auto&&... args)
    requires(is_field_v<Field_t> || is_field_tile_set_v<Field_t>)
{
    operate_field<Op>(dst, src, args...);
}

template<template<typename> typename Op, typename Field_t, typename PQ, std::size_t rank>
void operate(TensorField<Field_t, PQ, rank>& dst, TensorField<Field_t, PQ, rank> const& src,
             auto&&... args)
{
    for (std::size_t ci = 0; ci < dst.size(); ++ci)
        operate_field<Op>(dst[ci], src[ci], args...);
}


/** \brief calls fn(dst.data()[i], fields.data()[i]...) for every cell i of dst (and the
 * co-located fields, which must share dst's shape/tiling). Meant for callbacks that only
 * need the value at the same cell across several fields (no neighbor/projection access) -
 * unlike GridLayout::evalOnAnyBox, this works for tiled fields.
 */
template<typename Fn, typename Field_t, typename... Fields>
void operate_on_fields(Fn&& fn, Field_t& dst, Fields&... fields)
    requires(not is_field_tile_set_v<Field_t>)
{
    for (std::size_t i = 0; i < dst.size(); ++i)
        fn(dst.data()[i], fields.data()[i]...);
}

template<typename Fn, typename Field_t, typename... Fields>
void operate_on_fields(Fn&& fn, Field_t& dst, Fields&... fields)
    requires(is_field_tile_set_v<Field_t>)
{
    auto& dst_tiles = dst();
    for (std::size_t i = 0; i < dst_tiles.size(); ++i)
        operate_on_fields(fn, dst_tiles[i](), fields()[i]()...);
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ALGORITHM_HPP
