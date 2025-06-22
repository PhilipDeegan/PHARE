#ifndef PHARE_CORE_DEBUG_HPP
#define PHARE_CORE_DEBUG_HPP

#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"

#include <cassert>
#include <sstream>
#include <unordered_map>

// see std::unordered_map<...> debugging
#if !defined(PHARE_DEBUGGERINO)
#define PHARE_DEBUGGERINO 1
#endif


namespace PHARE::core
{

std::size_t constexpr inline static DEBUG_DIM = 2;

struct FieldActiveFunction
{
    double time;
    Box<int, DEBUG_DIM> box;
    std::string dst_name, src_name = dst_name;

    void operator()() const { PHARE_LOG_LINE_SS(""); }
};


std::unordered_map<std::string, FieldActiveFunction> static inline debugging{
    {"Simulator/advance/level/0/", {.005, Box<int, DEBUG_DIM>{{0, 0}, {1, 1}}, "EMAvg_E_x"}},
    {"Simulator/advance/level/0/SolverPPC/moveIons/0/",
     {.005, Box<int, DEBUG_DIM>{{-2, -2}, {2, 2}}, "HybridModel-HybridModel_sumVec_x"}},
    {"Simulator/advance/level/0/SolverPPC/moveIons/1/",
     {.005, Box<int, DEBUG_DIM>{{-2, -2}, {2, 2}}, "HybridModel-HybridModel_sumVec_x"}},
};


struct debug_scope;

struct Debuggerino
{
    static Debuggerino& INSTANCE();

    Debuggerino() = default;

    void settime(double const t) { time = t; }

    double time            = 0;
    debug_scope* stack_ptr = nullptr;
};




struct debug_scope
{
    debug_scope(std::string const& key);

    ~debug_scope();

    static debug_scope const* root_parent_from(auto const& self)
    {
        if (self.parent)
            return root_parent_from(*self.parent);
        return &self;
    }


    auto full_path() const
    {
        auto const self = this;
        auto iter       = root_parent_from(*this);
        std::stringstream ss;

        ss << iter->key;

        while (iter and iter != this)
            if ((iter = iter->child ? iter->child : nullptr))
                ss << iter->key;

        return ss.str();
    }

    std::string key;
    debug_scope* parent = nullptr;
    debug_scope* child  = nullptr;
};


void debug_fields(auto& dst, auto& src, auto& k, auto& v)
{
    auto& debugger           = Debuggerino::INSTANCE();
    auto& scope              = *debugger.stack_ptr;
    auto const path          = scope.full_path();
    auto const dst_ghost_box = shift(dst.amr_ghost_box, src.offset_ * -1);
    auto const src_ghost_box = src.amr_ghost_box;

    PHARE_LOG_LINE_SS("debugging: " << path << " / time: " << v.time
                                    << " / field: " << dst.field.name());

    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
    {
        auto const amr_idx = (*dst_it + dst_ghost_box.lower).as_signed();
        if (isIn(amr_idx, v.box))
        {
            auto const dv = to_string_with_precision(dst.field(*dst_it), 17);
            auto const sv = to_string_with_precision(src.field(*src_it), 17);

            std::cout << "amr_idx=" << amr_idx << " dst(" << dv << ") = src(" << sv << ")"
                      << std::endl;

            if (auto src_amr_idx = (*src_it + src_ghost_box.lower).as_signed();
                amr_idx != src_amr_idx)
            {
                PHARE_LOG_LINE_SS("amr_idx MISMATCH dst=" << amr_idx << " src=" << src_amr_idx);
            }
        }
    }
}

void debug_if_active_fields(auto const& dst, auto const& src)
{
    if constexpr (DEBUG_DIM == std::decay_t<decltype(dst)>::dimension)
    {
        assert(Debuggerino::INSTANCE().stack_ptr);
        auto& debugger           = Debuggerino::INSTANCE();
        debug_scope& scope       = *debugger.stack_ptr;
        auto const path          = scope.full_path();
        auto const dst_ghost_box = shift(dst.amr_ghost_box, src.offset_ * -1);
        auto const dst_amr_box   = shift(
            Box<int, DEBUG_DIM>{dst.lcl_box.lower.as_signed(), dst.lcl_box.upper.as_signed()},
            dst_ghost_box.lower);

        if (get_env("PHARE_PRINT_ALL_FIELDS"))
        {
            PHARE_LOG_LINE_SS(path << " / time: " << debugger.time
                                   << " / field: " << dst.field.name() << " " << dst_amr_box);
        }

        for (auto [k, v] : debugging)
            if (float_equals(debugger.time, v.time) and path == k and v.dst_name == dst.field.name()
                and v.box * dst_amr_box)
                debug_fields(dst, src, k, v);
    }
}



void debug_field(auto& field, auto& layout, auto& k, auto& v)
{
    auto& debugger           = Debuggerino::INSTANCE();
    auto& scope              = *debugger.stack_ptr;
    auto const path          = scope.full_path();
    auto const dst_ghost_box = layout.AMRGhostBoxFor(field.physicalQuantity());

    PHARE_LOG_LINE_SS("debugging: " << path << " / time: " << v.time
                                    << " / field: " << field.name());

    for (auto const& bix : layout.AMRToLocal(v.box))
    {
        auto const amr_idx = (bix + dst_ghost_box.lower).as_signed();
        if (isIn(amr_idx, v.box))
        {
            auto const dv = to_string_with_precision(field(*bix), 17);

            std::cout << "amr_idx=" << amr_idx << " dst(" << dv << ") " << std::endl;
        }
    }
}

void debug_if_active_field(auto const& field, auto const& layout)
{
    if constexpr (DEBUG_DIM == std::decay_t<decltype(field)>::dimension)
    {
        assert(Debuggerino::INSTANCE().stack_ptr);
        auto& debugger           = Debuggerino::INSTANCE();
        debug_scope& scope       = *debugger.stack_ptr;
        auto const path          = scope.full_path();
        auto const dst_ghost_box = layout.AMRGhostBoxFor(field.physicalQuantity());

        if (get_env("PHARE_PRINT_ALL_FIELDS"))
        {
            PHARE_LOG_LINE_SS(path << " / time: " << debugger.time << " / field: " << field.name()
                                   << " " << dst_ghost_box);
        }

        for (auto [k, v] : debugging)
            if (float_equals(debugger.time, v.time) and path == k and v.dst_name == field.name()
                and v.box * dst_ghost_box)
                debug_field(field, layout, k, v);
    }
}


void debug_all_fields(auto& views)
{
    using ResourcesManager_t = std::decay_t<decltype(*views.model().resourcesManager)>;
    using FieldData_t        = typename ResourcesManager_t::UserField_t::patch_data_type;
    auto const& rm           = *views.model().resourcesManager;

    for (auto [k, v] : rm.all_resources())
    {
        // PHARE_LOG_LINE_SS(k);

        for (auto& state : views)
        {
            auto data = state.patch->getPatchData(v.id);
            if (auto field_ptr = dynamic_cast<FieldData_t*>(data.get()))
                debug_if_active_field(*field_ptr->getPointer(), state.layout);
        }
    }
}


} // namespace PHARE::core

#if defined(PHARE_DEBUGGERINO) and PHARE_DEBUGGERINO

#define PHARE_DEBUG_SCOPE(key) PHARE::core::debug_scope _debug_scope_{key};
#define PHARE_DEBUG_TIME(time) PHARE::core::Debuggerino::INSTANCE().settime(time);
#define PHARE_DEBUG_FIELDS(...) PHARE::core::debug_if_active_fields(__VA_ARGS__);
#define PHARE_DEBUG_FIELD(...) PHARE::core::debug_if_active_field(__VA_ARGS__);
#define PHARE_DEBUG_ALL_FIELDS(...) PHARE::core::debug_all_fields(__VA_ARGS__);

#else // disabled

#define PHARE_DEBUG_SCOPE(key)
#define PHARE_DEBUG_TIME(time)
#define PHARE_DEBUG_FIELDS(...)
#define PHARE_DEBUG_FIELD(...)
#define PHARE_DEBUG_ALL_FIELDS(...)

#endif // PHARE_DEBUGGERINO

#endif /*PHARE_CORE_DEBUG_HPP*/
