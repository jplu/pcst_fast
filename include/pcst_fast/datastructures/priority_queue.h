#pragma once

#include <vector>
#include <utility>
#include <optional>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include "pcst_fast/pcst_types.h"

namespace cluster_approx {

template <typename ValueType, typename IndexType>
class PriorityQueue {
public:
    struct Element {
        ValueType value;
        IndexType index;
    };

    PriorityQueue() = default;

    void reserve(size_t max_elements, size_t max_index) {
        heap_.reserve(max_elements);
        pos_.assign(max_index, -1);
    }

    FORCE_INLINE [[nodiscard]] bool is_empty() const noexcept {
        return heap_.empty();
    }

    FORCE_INLINE [[nodiscard]] std::optional<std::pair<ValueType, IndexType>> get_min() const {
        if (heap_.empty()) return std::nullopt;
        return std::make_pair(heap_.front().value, heap_.front().index);
    }

    FORCE_INLINE std::optional<std::pair<ValueType, IndexType>> delete_min() {
        if (heap_.empty()) return std::nullopt;

        Element min_elem = heap_.front();
        pos_[min_elem.index] = -1;

        if (heap_.size() > 1) {
            heap_[0] = heap_.back();
            heap_.pop_back();
            pos_[heap_[0].index] = 0;
            sift_down(0);
        } else {
            heap_.pop_back();
        }
        return std::make_pair(min_elem.value, min_elem.index);
    }

    FORCE_INLINE void push_back_fast(ValueType value, IndexType index) {
        size_t idx = static_cast<size_t>(index);
        pos_[idx] = static_cast<int32_t>(heap_.size());
        heap_.push_back({value, index});
    }

    FORCE_INLINE void build_heap() {
        int32_t n = static_cast<int32_t>(heap_.size());
        if (n <= 1) return;
        for (int32_t i = (n - 2) / 4; i >= 0; --i) {
            sift_down(i);
        }
    }

    FORCE_INLINE void insert_or_update(ValueType value, IndexType index) {
        size_t idx = static_cast<size_t>(index);
        assert(idx < pos_.size());
        int32_t p = pos_[idx];

        if (p == -1) {
            p = static_cast<int32_t>(heap_.size());
            heap_.push_back({value, index});
            pos_[idx] = p;
            sift_up(p);
        } else {
            if (value < heap_[p].value) {
                heap_[p].value = value;
                sift_up(p);
            } else if (value > heap_[p].value) {
                heap_[p].value = value;
                sift_down(p);
            }
        }
    }

    FORCE_INLINE void insert(ValueType value, IndexType index) {
        insert_or_update(value, index);
    }

    FORCE_INLINE void decrease_key(ValueType new_value, IndexType index) {
        insert_or_update(new_value, index);
    }

    FORCE_INLINE void delete_element(IndexType index) {
        size_t idx = static_cast<size_t>(index);
        assert(idx < pos_.size());
        int32_t p = pos_[idx];
        if (p == -1) return;

        int32_t last_idx = static_cast<int32_t>(heap_.size()) - 1;
        if (p == last_idx) {
            pos_[idx] = -1;
            heap_.pop_back();
        } else {
            heap_[p] = heap_.back();
            heap_.pop_back();
            pos_[heap_[p].index] = p;
            pos_[idx] = -1;

            if (p > 0 && heap_[p].value < heap_[(p - 1) / 4].value) {
                sift_up(p);
            } else {
                sift_down(p);
            }
        }
    }

private:
    std::vector<Element> heap_;
    std::vector<int32_t> pos_;

    FORCE_INLINE void sift_up(int32_t i) {
        if (i == 0) return;
        Element val = heap_[i];
        while (i > 0) {
            int32_t p = (i - 1) / 4;
            if (val.value < heap_[p].value) {
                heap_[i] = heap_[p];
                pos_[heap_[i].index] = i;
                i = p;
            } else {
                break;
            }
        }
        heap_[i] = val;
        pos_[val.index] = i;
    }

    FORCE_INLINE void sift_down(int32_t i) {
        int32_t n = static_cast<int32_t>(heap_.size());
        if (n <= 1) return;
        Element val = heap_[i];

        while (true) {
            int32_t first_child = 4 * i + 1;
            if (first_child >= n) break;

            int32_t min_child = first_child;
            ValueType min_val = heap_[first_child].value;

            int32_t limit = std::min(first_child + 4, n);
            for (int32_t child = first_child + 1; child < limit; ++child) {
                if (heap_[child].value < min_val) {
                    min_val = heap_[child].value;
                    min_child = child;
                }
            }

            if (min_val < val.value) {
                heap_[i] = heap_[min_child];
                pos_[heap_[i].index] = i;
                i = min_child;
            } else {
                break;
            }
        }
        heap_[i] = val;
        pos_[val.index] = i;
    }
};

}