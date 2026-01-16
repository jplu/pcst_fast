#pragma once

#include <vector>
#include <utility>
#include <optional>
#include <cassert>
#include <algorithm>

namespace cluster_approx {

/**
 * @brief A binary heap priority queue optimized with an index map for fast decrease-key.
 */
template <typename ValueType, typename IndexType>
class PriorityQueue {
  public:
    struct Element {
        ValueType value;
        IndexType index;
    };

    PriorityQueue() = default;

    [[nodiscard]] bool is_empty() const noexcept {
        return heap_.empty();
    }

    [[nodiscard]] std::optional<std::pair<ValueType, IndexType>> get_min() const {
        if (heap_.empty()) return std::nullopt;
        return std::make_pair(heap_.front().value, heap_.front().index);
    }

    [[nodiscard]] std::optional<std::pair<ValueType, IndexType>> delete_min() {
        if (heap_.empty()) return std::nullopt;

        Element min_elem = heap_.front();
        remove_node(0);
        return std::make_pair(min_elem.value, min_elem.index);
    }

    void insert(ValueType value, IndexType index) {
        assert(index >= 0);
        if (static_cast<size_t>(index) >= position_map_.size()) {
            position_map_.resize(static_cast<size_t>(index) + 1, -1);
        }

        if (position_map_[index] != -1) {
            remove_node(position_map_[index]);
        }

        heap_.push_back({value, index});
        int i = static_cast<int>(heap_.size()) - 1;
        position_map_[index] = i;
        sift_up(i);
    }

    void decrease_key(ValueType new_value, IndexType index) {
        assert(static_cast<size_t>(index) < position_map_.size() && position_map_[index] != -1);
        int i = position_map_[index];
        heap_[i].value = new_value;
        sift_up(i);
    }

    void delete_element(IndexType index) {
        if (static_cast<size_t>(index) < position_map_.size()) {
            int i = position_map_[index];
            if (i != -1) {
                remove_node(i);
            }
        }
    }

  private:
    std::vector<Element> heap_;
    std::vector<int> position_map_; // Maps external IndexType to heap position

    void sift_up(int i) {
        while (i > 0) {
            int p = (i - 1) / 2;
            if (heap_[i].value < heap_[p].value) {
                swap_nodes(i, p);
                i = p;
            } else {
                break;
            }
        }
    }

    void sift_down(int i) {
        int n = static_cast<int>(heap_.size());
        while (true) {
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            int smallest = i;

            if (left < n && heap_[left].value < heap_[smallest].value)
                smallest = left;
            if (right < n && heap_[right].value < heap_[smallest].value)
                smallest = right;

            if (smallest != i) {
                swap_nodes(i, smallest);
                i = smallest;
            } else {
                break;
            }
        }
    }

    void swap_nodes(int i, int j) {
        std::swap(heap_[i], heap_[j]);
        position_map_[heap_[i].index] = i;
        position_map_[heap_[j].index] = j;
    }

    void remove_node(int i) {
        int last = static_cast<int>(heap_.size()) - 1;
        if (i != last) {
            swap_nodes(i, last);
            position_map_[heap_[last].index] = -1;
            heap_.pop_back();
            sift_up(i);   
            sift_down(i); 
        } else {
            position_map_[heap_[last].index] = -1;
            heap_.pop_back();
        }
    }
};

}