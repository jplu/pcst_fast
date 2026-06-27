#pragma once

#include <vector>
#include <utility>
#include <optional>
#include <cassert>
#include <algorithm>
#include <cstdint>

namespace cluster_approx {

/**
 * @brief A 4-ary heap priority queue using lazy deletions.
 *
 * This design eliminates the expensive lookup mapping and active decrease-key operations.
 * Keys updates are handled as duplicate pushes with incremented version IDs.
 * Stale duplicates are cleared from the top of the heap only when querying or popping.
 * A 4-ary layout matches CPU cache lines better than traditional binary heaps.
 */
template <typename ValueType, typename IndexType>
class PriorityQueue {
public:
    struct Element {
        ValueType value;
        IndexType index;
        uint32_t version;

        bool operator<(const Element& other) const noexcept {
            return value < other.value;
        }
        bool operator>(const Element& other) const noexcept {
            return value > other.value;
        }
    };

    PriorityQueue() = default;

    [[nodiscard]] bool is_empty() const noexcept {
        purge_stale();
        return heap_.empty();
    }

    [[nodiscard]] std::optional<std::pair<ValueType, IndexType>> get_min() const {
        purge_stale();
        if (heap_.empty()) return std::nullopt;
        return std::make_pair(heap_.front().value, heap_.front().index);
    }

    [[nodiscard]] std::optional<std::pair<ValueType, IndexType>> delete_min() {
        purge_stale();
        if (heap_.empty()) return std::nullopt;

        Element min_elem = heap_.front();
        pop_front();
        return std::make_pair(min_elem.value, min_elem.index);
    }

    void insert(ValueType value, IndexType index) {
        assert(index >= 0);
        size_t idx = static_cast<size_t>(index);
        if (idx >= versions_.size()) {
            versions_.resize(idx + 1, 0);
        }

        // Increment the current version to invalidate prior heap instances of this index
        uint32_t new_version = ++versions_[idx];
        heap_.push_back({value, index, new_version});
        sift_up(heap_.size() - 1);
    }

    void decrease_key(ValueType new_value, IndexType index) {
        // Under lazy deletion, updating a key is treated as inserting with a fresh version
        insert(new_value, index);
    }

    void delete_element(IndexType index) {
        size_t idx = static_cast<size_t>(index);
        if (idx < versions_.size()) {
            // Incrementing the target index's version automatically invalidates its active heap nodes
            versions_[idx]++;
        }
    }

private:
    mutable std::vector<Element> heap_;
    mutable std::vector<uint32_t> versions_;

    void purge_stale() const {
        while (!heap_.empty()) {
            const auto& top = heap_.front();
            size_t idx = static_cast<size_t>(top.index);
            if (idx < versions_.size() && top.version == versions_[idx]) {
                break; // Found valid active event
            }
            pop_front(); // Evict stale event
        }
    }

    void pop_front() const {
        if (heap_.empty()) return;
        if (heap_.size() > 1) {
            heap_[0] = std::move(heap_.back());
            heap_.pop_back();
            sift_down(0);
        } else {
            heap_.pop_back();
        }
    }

    void sift_up(size_t i) const {
        while (i > 0) {
            size_t p = (i - 1) / 4;
            if (heap_[i] < heap_[p]) {
                std::swap(heap_[i], heap_[p]);
                i = p;
            } else {
                break;
            }
        }
    }

    void sift_down(size_t i) const {
        size_t n = heap_.size();
        while (true) {
            size_t first_child = 4 * i + 1;
            if (first_child >= n) break;

            size_t smallest = i;
            // Scan through up to 4 children of the current node
            for (size_t c = 0; c < 4; ++c) {
                size_t child = first_child + c;
                if (child < n) {
                    if (heap_[child] < heap_[smallest]) {
                        smallest = child;
                    }
                } else {
                    break;
                }
            }

            if (smallest != i) {
                std::swap(heap_[i], heap_[smallest]);
                i = smallest;
            } else {
                break;
            }
        }
    }
};

}