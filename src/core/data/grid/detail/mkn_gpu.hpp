
#ifndef PHARE_CORE_GRID_DETAIL_MKN_GPU_H
#define PHARE_CORE_GRID_DETAIL_MKN_GPU_H

namespace PHARE::core::mkn_gpu
{
struct GridLayout
{
    template<typename Fn, typename Index, typename... Args>
    static void evalOnBox(Fn& fn, std::vector<Index>& indexes, Args&... args)
    {
        // MKN_GPU::resources::Cuda res;
        // auto tuple     = std::make_tuple(args...); // copy for copy :(
        // auto gpu_tuple = allocate_copy(res, tuple);
        // auto d_indexes = allocate_copy(res, indexes);
        // exec(
        //     [=] MKN_GPU_DEVICE(int i) mutable {
        //         std::apply(
        //             [=](auto&... targs) mutable { //
        //                 fn(d_indexes[i], targs...);
        //             },
        //             *gpu_tuple);
        //     },
        //     res, indexes.size());
        // deallocate(res, gpu_tuple, d_indexes);
    }
};


} // namespace PHARE::core::mkn_gpu
#endif /* PHARE_CORE_GRID_DETAIL_MKN_GPU_H */
