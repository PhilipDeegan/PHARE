#ifndef PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP
#define PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP

#include "core/def/phare_config.hpp"

#if PHARE_HAVE_MKN_GPU
// #include "mkn_boris.hpp"
// #include "mkn_boris_capture.hpp"
#include "tile_boris_capture.hpp"
#endif // PHARE_HAVE_MKN_GPU

#include <stdexcept>

namespace PHARE::core::boris
{


template<typename... Args>
auto mover(Args&&... args)
{
    if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(return MultiBoris{args...});
    }
    else
        throw std::runtime_error("Vector::copy NO ALTERNATIVE");
}


} // namespace PHARE::core::boris


#endif /*PHARE_CORE_PUSHER_BORIS_BORIS_MOVE_HPP*/
