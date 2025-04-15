#pragma once

#include <set>
#include <vector>
#include <utility>
#include <optional>
#include <limits>
#include <cassert>

namespace cluster_approx {

/**
 * @brief A priority queue implementation supporting decrease-key and delete-element operations.
 *
 * Internally uses std::set for ordering and std::vector to store iterators
 * for efficient access to elements by index for updates/deletions. Assumes
 * IndexType is a non-negative integer used for indexing the iterator vector.
 *
 * @tparam ValueType The type of the priority values (keys). Must be comparable.
 * @tparam IndexType The type of the index associated with each value (payload). Must be suitable as an index for std::vector.
 */
template <typename ValueType, typename IndexType>
class PriorityQueue {
  private:
    using SetPair = std::pair<ValueType, IndexType>;
    using SetIterator = typename std::set<SetPair>::iterator;

  public:
    PriorityQueue() = default;

    /**
     * @brief Checks if the priority queue is empty.
     * @return True if the queue contains no elements, false otherwise.
     */
    [[nodiscard]] bool is_empty() const noexcept {
        return sorted_set_.empty();
    }

    /**
     * @brief Gets the minimum element (value and index) without removing it.
     * @return An std::optional containing the pair (value, index) of the minimum element if the queue is not empty, otherwise std::nullopt.
     */
    [[nodiscard]] std::optional<SetPair> get_min() const {
        if (sorted_set_.empty()) {
            return std::nullopt;
        }

        return *sorted_set_.begin();
    }

    /**
     * @brief Removes and returns the minimum element (value and index).
     * @return An std::optional containing the pair (value, index) of the minimum element if the queue was not empty, otherwise std::nullopt.
     */
    [[nodiscard]] std::optional<SetPair> delete_min() {
        if (sorted_set_.empty()) {
            return std::nullopt;
        }
        SetPair min_pair = *sorted_set_.begin();
        sorted_set_.erase(sorted_set_.begin());

        return min_pair;
    }

    /**
     * @brief Inserts a new element or updates an existing element with the given index.
     * If an element with the same index already exists, it's effectively updated (removed and re-inserted).
     * @param value The priority value of the element.
     * @param index The index associated with the element. Must be non-negative.
     */
    void insert(ValueType value, IndexType index) {
        assert(index >= 0 && "Index must be non-negative.");
        if (static_cast<size_t>(index) >= index_to_iterator_.size()) {

            index_to_iterator_.resize(static_cast<size_t>(index) + 1, sorted_set_.end());
        } else {

            SetIterator existing_it = index_to_iterator_[index];
            if (existing_it != sorted_set_.end()) {

                sorted_set_.erase(existing_it);
            }
        }

        auto insert_result = sorted_set_.insert({value, index});
        assert(insert_result.second && "Insertion into set failed unexpectedly.");
        index_to_iterator_[index] = insert_result.first;
    }

    /**
     * @brief Decreases the priority value of an existing element.
     * If the new value is not strictly smaller, the behavior might be unexpected
     * depending on the std::set implementation details (element might not move).
     * The original code replaces the element regardless. We mimic this.
     * @param new_value The new, smaller priority value.
     * @param index The index of the element to update. Must exist in the queue.
     */
    void decrease_key(ValueType new_value, IndexType index) {
        assert(static_cast<size_t>(index) < index_to_iterator_.size() && "Index out of bounds.");
        SetIterator it = index_to_iterator_[index];
        assert(it != sorted_set_.end() && "Attempting to decrease key for non-existent or deleted index.");

        sorted_set_.erase(it);
        auto insert_result = sorted_set_.insert({new_value, index});
        assert(insert_result.second && "Re-insertion failed in decrease_key.");
        index_to_iterator_[index] = insert_result.first;
    }

    /**
     * @brief Removes an element from the priority queue by its index.
     * Does nothing if the index is out of bounds or the element is not currently in the queue.
     * @param index The index of the element to remove.
     */
    void delete_element(IndexType index) {
        assert(index >= 0 && "Index must be non-negative.");
        if (static_cast<size_t>(index) < index_to_iterator_.size()) {
            SetIterator it = index_to_iterator_[index];
            if (it != sorted_set_.end()) {
                sorted_set_.erase(it);

                index_to_iterator_[index] = sorted_set_.end();
            }
        }

    }

  private:
    std::set<SetPair> sorted_set_;

    std::vector<SetIterator> index_to_iterator_;
};

}