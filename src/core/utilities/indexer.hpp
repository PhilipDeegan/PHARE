#ifndef PHARE_INDEXER_H
#define PHARE_INDEXER_H

#include "core/data/particles/particle_array_def.hpp"
#include "core/def.hpp"
#include "core/utilities/span.hpp"
// #include "core/utilities/types.hpp"

#include <array>
#include <vector>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <unordered_map>

namespace PHARE::core
{

template<typename Indexes>
struct IndexerBase
{
    IndexerBase() = default;
    IndexerBase(auto const& data)
        : indexes_{data}
    {
    }

    NO_DISCARD auto data() const _PHARE_ALL_FN_ { return indexes_.data(); }
    NO_DISCARD auto data() _PHARE_ALL_FN_ { return indexes_.data(); }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return indexes_.size(); }


    NO_DISCARD bool is_indexed(std::size_t itemIndex)
    {
        return std::end(indexes_) != std::find(std::begin(indexes_), std::end(indexes_), itemIndex);
    }


    Indexes indexes_;
};


template<auto storage_mode>
struct IndexerStorage;


template<>
struct IndexerStorage<StorageMode::SPAN> : public IndexerBase<Span<std::size_t>>
{
    using Super = IndexerBase<Span<std::size_t>>;
    using Super::indexes_;

    IndexerStorage() = default;
    // IndexerStorage(auto const& span)
    //     : Super{span}
    // {
    // }

    void set_from(auto& indexer)
    {
        indexes_.ptr = indexer.data(); // = Span<std::size_t>{indexer.data(), indexer.size()};
        indexes_.s   = indexer.size();
    }
};

template<>
struct IndexerStorage<StorageMode::VECTOR> : public IndexerBase<std::vector<std::size_t>>
{
    using Super = IndexerBase<std::vector<std::size_t>>;
    using Super::indexes_;

    IndexerStorage() = default;
    IndexerStorage(auto const& data)
        : Super{data}
    {
    }

    void add(std::size_t itemIndex) { indexes_.push_back(itemIndex); }
    void remove(std::size_t itemIndex)
    {
        auto it = std::find(std::begin(indexes_), std::end(indexes_), itemIndex);
        if (it != std::end(indexes_))
        {
            indexes_.erase(it);
        }
        assert(!is_indexed(itemIndex));
    }

    // empty the bucketlist, but leaves the capacity untouched
    void empty() { indexes_.resize(0); };
    NO_DISCARD std::size_t capacity() const { return indexes_.capacity(); }
    void resize(std::size_t const s) { indexes_.resize(s); };
    void clear() { indexes_.clear(); }

    // std::vector<std::size_t> indexes_;
};


template<auto storage_mode = StorageMode::VECTOR>
class Indexer : public IndexerStorage<storage_mode>
{
    using Super = IndexerStorage<storage_mode>;
    using Super::indexes_;

public:
    using Super::data;
    using Super::size;

    Indexer()                                = default;
    Indexer(Indexer const& other)            = default;
    Indexer(Indexer&& other)                 = default;
    Indexer& operator=(Indexer const& other) = default;
    Indexer& operator=(Indexer&& other)      = default;




    // reallocate bucketlist memory if more empty space than max_empty
    // void trim(std::size_t max_empty);

    // to use if an item in an indexed array is moved at another index
    void updateIndex(std::size_t oldIndex, std::size_t newIndex)
    {
        //
        auto it = std::find(std::begin(indexes_), std::end(indexes_), oldIndex);
        if (it != std::end(indexes_))
        {
            *it = newIndex;
        }
    }



    NO_DISCARD bool is_empty() const { return indexes_.size() == 0; }

    NO_DISCARD std::size_t size() const { return indexes_.size(); }
    NO_DISCARD auto begin() { return indexes_.begin(); }
    NO_DISCARD auto begin() const { return indexes_.begin(); }
    NO_DISCARD auto cbegin() const { return indexes_.begin(); }
    NO_DISCARD auto end() { return indexes_.end(); }
    NO_DISCARD auto end() const { return indexes_.end(); }
    NO_DISCARD auto cend() const { return indexes_.end(); }
    void sort() { std::sort(indexes_.begin(), indexes_.end()); }

    NO_DISCARD auto& operator[](std::size_t const idx) { return indexes_[idx]; }
    NO_DISCARD auto& operator[](std::size_t const idx) const { return indexes_[idx]; }

    NO_DISCARD auto& back() { return indexes_.back(); }
    NO_DISCARD auto& front() const { return indexes_.front(); }
    void swap(auto const& a, auto const& b) { std::swap(indexes_[a], indexes_[b]); }
};


// ==========================================================




} // namespace PHARE::core


#endif
