#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <functional>
#include <type_traits>
#include <cstddef>
#include <algorithm>
#include <concepts>
#include <stdexcept> // For potential exceptions
#include <optional>

namespace cluster_approx {
    template<typename T>
    concept Indexable = std::is_integral_v<T> || std::is_enum_v<T>;

    namespace internal {
        // Helper to get underlying type for enums, or the type itself otherwise
        template <typename T, bool = std::is_enum_v<T>>
        struct get_underlying_or_self_impl { using type = T; };

        template <typename T>
        struct get_underlying_or_self_impl<T, true> { using type = std::underlying_type_t<T>; };

        template <Indexable T>
        using get_underlying_or_self_t = typename get_underlying_or_self_impl<T>::type;

        // Comparator for heap elements (pair<Value, Index>)
        // Breaks ties using index to ensure deterministic behavior.
        template <typename V, Indexable I, typename C>
            requires std::strict_weak_order<C, V, V>
        struct HeapElementCompare {
            [[no_unique_address]] C value_comp; // Use attribute for potential space saving
            // std::less<> for index comparison ensures total order if values are equal
            std::less<I> index_comp {};

            explicit HeapElementCompare(const C& comp = C()) noexcept(std::is_nothrow_copy_constructible_v<C>)
                : value_comp(comp) {}

            [[nodiscard]] constexpr bool operator()(const std::pair<V, I>& lhs, const std::pair<V, I>& rhs) const
                noexcept(noexcept(value_comp(lhs.first, rhs.first)) && noexcept(index_comp(lhs.second, rhs.second)))
            {
                // Primary comparison based on value
                if (value_comp(lhs.first, rhs.first)) return true;
                if (value_comp(rhs.first, lhs.first)) return false;
                // Secondary comparison (tie-breaker) based on index
                return index_comp(lhs.second, rhs.second);
            }
        };
    } // namespace internal

    /**
     * @brief A binary heap based Priority Queue supporting decrease_key and delete_element by index.
     *
     * Maintains a mapping from IndexType to the element's position in the heap.
     *
     * @tparam ValueType The type of the values stored (priorities). Must be totally_ordered and movable.
     * @tparam IndexType The type used to index elements (e.g., int, size_t, enum). Must be Indexable and movable.
     * @tparam Compare   The comparison function object type (defaults to std::less for a min-heap).
     *
     * @note This class is NOT thread-safe. Concurrent access to the same instance
     *       requires external synchronization (e.g., std::mutex). However, multiple
     *       instances can be used independently in different threads.
     */
    template <typename ValueType,
              Indexable IndexType,
              typename Compare = std::less<ValueType>>
        requires std::strict_weak_order<Compare, ValueType, ValueType>
              && std::totally_ordered<ValueType> // Required by some operations implicitly
              && std::movable<ValueType>
              && std::movable<IndexType>
    class PriorityQueue {
    private:
        using PositionType = std::size_t;
        using UnderlyingIndexType = internal::get_underlying_or_self_t<IndexType>;
        using HeapElement = std::pair<ValueType, IndexType>;
        using Comp = internal::HeapElementCompare<ValueType, IndexType, Compare>;

        // Constants
        static constexpr PositionType INVALID_POS = std::numeric_limits<PositionType>::max();
        static constexpr PositionType ROOT_POS = 0;

        // Member Variables
        std::vector<HeapElement> heap_data_; // Stores {value, index} pairs in heap order
        std::vector<PositionType> index_to_heap_pos_; // Maps IndexType -> PositionType in heap_data_
        [[no_unique_address]] Comp heap_comp_; // The comparator instance

        // --- Heap Structure Helpers ---

        [[nodiscard]] static constexpr PositionType parent(PositionType i) noexcept { return (i - 1) / 2; }
        [[nodiscard]] static constexpr PositionType left(PositionType i) noexcept { return (i * 2) + 1; }
        [[nodiscard]] static constexpr PositionType right(PositionType i) noexcept { return (i * 2) + 2; }

        // --- Index Mapping Helpers ---

        [[nodiscard]] static constexpr PositionType to_pos_index(IndexType index) noexcept {
            // Cast IndexType (potentially enum or signed) to UnderlyingIndexType (unsigned or signed integral)
            // then to PositionType (size_t). Assumes non-negative indices if Underlying is signed.
             // A check for negative values might be useful if IndexType is signed.
             // static_assert(std::is_unsigned_v<UnderlyingIndexType> || std::is_signed_v<UnderlyingIndexType>, "..."); // Already checked in original
             UnderlyingIndexType underlying_idx = static_cast<UnderlyingIndexType>(index);
             // Add assertion for signed types if negative indices are invalid
             // assert( !std::is_signed_v<UnderlyingIndexType> || underlying_idx >= 0);
            return static_cast<PositionType>(underlying_idx);
        }

        // Resizes index_to_heap_pos_ if needed to accommodate index.
        void ensure_index_capacity(IndexType index) {
            const PositionType pos_index = to_pos_index(index);
            if (pos_index >= index_to_heap_pos_.size()) {
                 // Grow geometrically, ensuring enough space for pos_index
                PositionType current_size = index_to_heap_pos_.size();
                PositionType desired_size = std::max(current_size + (current_size >> 1) + 1, pos_index + 1);
                index_to_heap_pos_.resize(desired_size, INVALID_POS);
            }
        }

        // Sets the heap position for a given index (unsafe: assumes index is valid).
        constexpr void set_pos_unsafe(IndexType index, PositionType heap_pos) noexcept {
             // Assumes ensure_index_capacity has been called or index is within bounds
            index_to_heap_pos_[to_pos_index(index)] = heap_pos;
        }

        // Gets the heap position for a given index (unsafe: assumes index is valid).
        [[nodiscard]] constexpr PositionType get_pos_unsafe(IndexType index) const noexcept {
             // Assumes ensure_index_capacity has been called or index is within bounds
            return index_to_heap_pos_[to_pos_index(index)];
        }

        // Checks if a heap position `pos` is valid and consistent with the index map at `pos_index`.
         [[nodiscard]] constexpr bool is_valid_and_consistent(PositionType pos, PositionType pos_index) const noexcept {
            return pos != INVALID_POS                   // Position must be valid
                && pos < heap_data_.size()            // Position must be within heap bounds
                && pos_index < index_to_heap_pos_.size() // Index must be within map bounds
                && to_pos_index(heap_data_[pos].second) == pos_index; // Heap element index must match map index
        }

        // Checks if index `index` currently exists in the heap.
        [[nodiscard]] bool contains(IndexType index) const noexcept {
            const PositionType pos_index = to_pos_index(index);
            if (pos_index >= index_to_heap_pos_.size()) {
                return false;
            }
            const PositionType heap_pos = index_to_heap_pos_[pos_index];
            return is_valid_and_consistent(heap_pos, pos_index);
        }


        // --- Core Heap Operations (Private) ---

        // Moves element at index `i` up the heap to restore heap property.
        void sift_up(PositionType i) noexcept(
            noexcept(std::move(heap_data_[i])) &&
            noexcept(heap_comp_(std::declval<HeapElement&>(), std::declval<HeapElement&>())) &&
            noexcept(set_pos_unsafe(std::declval<IndexType&>(), std::declval<PositionType&>()))
        )
        {
            if (i == ROOT_POS || i >= heap_data_.size()) return;

            HeapElement target = std::move(heap_data_[i]); // Move the element out
            PositionType current_pos = i;

            // Loop while not at root and target is smaller than parent
            while (current_pos > ROOT_POS) {
                PositionType p_pos = parent(current_pos);
                if (!heap_comp_(target, heap_data_[p_pos])) {
                    break; // Found correct position
                }

                // Move parent down
                heap_data_[current_pos] = std::move(heap_data_[p_pos]);
                set_pos_unsafe(heap_data_[current_pos].second, current_pos); // Update moved parent's index

                current_pos = p_pos; // Move up
            }

            // Place target in its final position
            heap_data_[current_pos] = std::move(target);
            set_pos_unsafe(heap_data_[current_pos].second, current_pos);
        }

        // Moves element at index `i` down the heap to restore heap property.
        void sift_down(PositionType i) noexcept(
            noexcept(std::move(heap_data_[i])) &&
            noexcept(heap_comp_(std::declval<HeapElement&>(), std::declval<HeapElement&>())) &&
            noexcept(set_pos_unsafe(std::declval<IndexType&>(), std::declval<PositionType&>()))
        )
        {
            const PositionType heap_size = heap_data_.size();
            if (heap_size <= 1 || i >= heap_size) return;

            HeapElement target = std::move(heap_data_[i]); // Move the element out
            PositionType current_pos = i;

            // Loop while node has at least a left child
            while (true) {
                PositionType child_pos = left(current_pos);
                if (child_pos >= heap_size) {
                    break; // No children, stop
                }

                PositionType smaller_child_pos = find_smaller_child(current_pos, heap_size);

                // If the smaller child is smaller than the target, swap and continue down
                if (heap_comp_(heap_data_[smaller_child_pos], target)) {
                    heap_data_[current_pos] = std::move(heap_data_[smaller_child_pos]);
                    set_pos_unsafe(heap_data_[current_pos].second, current_pos);
                    current_pos = smaller_child_pos; // Move down
                } else {
                    break; // Target is smaller or equal, stop
                }
            }

            // Place target in its final position
            heap_data_[current_pos] = std::move(target);
            set_pos_unsafe(heap_data_[current_pos].second, current_pos);
        }

        // Helper for sift_down: Finds the index of the smaller child of node i.
        [[nodiscard]] PositionType find_smaller_child(PositionType i, PositionType heap_size) const noexcept {
             PositionType l_pos = left(i);
             PositionType r_pos = right(i);
             PositionType smaller_child_pos = l_pos; // Assume left is smaller initially

             // If right child exists and is smaller than left child
             if (r_pos < heap_size && heap_comp_(heap_data_[r_pos], heap_data_[l_pos])) {
                 smaller_child_pos = r_pos;
             }
             return smaller_child_pos;
        }


        // Removes the element at heap position `pos`, replaces it with the last element,
        // and sifts appropriately. Assumes `pos` is valid.
        void remove_at(PositionType pos) noexcept (
             noexcept(std::move(heap_data_[pos])) &&
             noexcept(sift_up(pos)) &&
             noexcept(sift_down(pos))
             )
         {
            const PositionType pos_index = to_pos_index(heap_data_[pos].second);
            index_to_heap_pos_[pos_index] = INVALID_POS; // Mark as removed in map

            const PositionType last_pos = heap_data_.size() - 1;

            if (pos == last_pos) {
                heap_data_.pop_back(); // Removing the last element is easy
            } else {
                 // Move the last element to the position being removed
                heap_data_[pos] = std::move(heap_data_[last_pos]);
                heap_data_.pop_back();
                set_pos_unsafe(heap_data_[pos].second, pos); // Update moved element's index

                // Restore heap property: try sifting up first, then down if needed
                // (only one of them will actually do work)
                if (pos > ROOT_POS && heap_comp_(heap_data_[pos], heap_data_[parent(pos)])) {
                    sift_up(pos);
                } else {
                    sift_down(pos); // Will also handle pos == ROOT_POS correctly
                }
            }
        }


    public:
        // --- Constructors, Destructor, Assignment ---

        explicit PriorityQueue(const Compare& compare = Compare()) noexcept(noexcept(Comp(compare)))
            : heap_comp_(compare) {}

        // Defaulted move operations are fine as members handle moves correctly.
        PriorityQueue(PriorityQueue&&) noexcept = default;
        PriorityQueue& operator=(PriorityQueue&&) noexcept = default;

        // Delete copy operations as they are complex and potentially expensive.
        PriorityQueue(const PriorityQueue&) = delete;
        PriorityQueue& operator=(const PriorityQueue&) = delete;

        ~PriorityQueue() = default; // Default destructor is sufficient

        // --- Public Interface ---

        [[nodiscard]] bool is_empty() const noexcept {
            return heap_data_.empty();
        }

        [[nodiscard]] std::size_t size() const noexcept {
            return heap_data_.size();
        }

        /**
         * @brief Gets the minimum element (value and index) without removing it.
         * @param value Pointer to store the minimum value.
         * @param index Pointer to store the index of the minimum element.
         * @return True if the heap is not empty, false otherwise.
         */
        [[nodiscard]] bool get_min(ValueType* value, IndexType* index) const
            noexcept(std::is_nothrow_copy_assignable_v<ValueType> && std::is_nothrow_copy_assignable_v<IndexType>)
        {
            if (is_empty()) {
                // Assign default values only if pointers are non-null
                if (value) *value = ValueType{};
                if (index) *index = IndexType{};
                return false;
            }
            // Check output pointers validity
            if (!value || !index) {
                // Or throw std::invalid_argument("Output pointers cannot be null");
                 return false; // Keep original behavior
            }

            *value = heap_data_[ROOT_POS].first;
            *index = heap_data_[ROOT_POS].second;
            return true;
        }

        /**
         * @brief Gets the minimum element (value and index) without removing it.
         * @return An optional containing the pair {value, index} if the heap is not empty,
         *         std::nullopt otherwise.
         */
        [[nodiscard]] std::optional<HeapElement> get_min() const {
             if(is_empty()) {
                 return std::nullopt;
             }
             return heap_data_[ROOT_POS];
        }


        /**
         * @brief Removes the minimum element from the heap.
         * @param value Pointer to store the value of the removed element.
         * @param index Pointer to store the index of the removed element.
         * @return True if an element was removed, false if the heap was empty.
         */
        [[nodiscard]] bool delete_min(ValueType* value, IndexType* index)
            noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_down(ROOT_POS)))
        {
            if (is_empty()) {
                if (value) *value = ValueType{};
                if (index) *index = IndexType{};
                return false;
            }
             if (!value || !index) {
                 // Or throw std::invalid_argument("Output pointers cannot be null");
                 return false; // Keep original behavior
             }

            // Retrieve min element data
            *value = std::move(heap_data_[ROOT_POS].first); // Move value out
            *index = heap_data_[ROOT_POS].second;

            // Remove element at root (handles swap, pop, sift)
            remove_at(ROOT_POS);

            return true;
        }

       /**
        * @brief Removes the minimum element from the heap.
        * @return An optional containing the {value, index} of the removed element,
        *         or std::nullopt if the heap was empty.
        */
        std::optional<HeapElement> delete_min()
             noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_down(ROOT_POS)))
        {
            if (is_empty()) {
                 return std::nullopt;
            }
            HeapElement min_element = std::move(heap_data_[ROOT_POS]); // Move element out
            remove_at(ROOT_POS);
            return min_element;
        }


        /**
         * @brief Inserts a new element or updates the value if the index already exists.
         * If the index exists, its value is updated (can be increase or decrease),
         * and the heap property is restored. If it doesn't exist, it's inserted.
         *
         * @param value The value of the element.
         * @param index The index of the element.
         */
        void insert_or_update(ValueType value, IndexType index)
             noexcept(noexcept(ensure_index_capacity(index)) &&
                      noexcept(delete_element(index)) && // Used internally by original insert
                      noexcept(heap_data_.emplace_back(std::move(value), index)) &&
                      noexcept(sift_up(0)) && noexcept(sift_down(0))) // For update case
        {
            ensure_index_capacity(index);
            const PositionType pos_index = to_pos_index(index);
            const PositionType current_heap_pos = index_to_heap_pos_[pos_index]; // Safe after ensure_capacity

            if (is_valid_and_consistent(current_heap_pos, pos_index)) {
                 // --- Update existing element ---
                 ValueType& current_value = heap_data_[current_heap_pos].first;
                 bool value_changed = !(current_value == value); // Requires operator== for ValueType
                 bool decreased = value_changed && heap_comp_.value_comp(value, current_value);
                 bool increased = value_changed && heap_comp_.value_comp(current_value, value);

                 if (value_changed) {
                     heap_data_[current_heap_pos].first = std::move(value);
                     if (decreased) {
                         sift_up(current_heap_pos);
                     } else if (increased) {
                         sift_down(current_heap_pos);
                     }
                 } // else: value is the same, do nothing
            } else {
                 // --- Insert new element ---
                 if (pos_index < index_to_heap_pos_.size()) { // Ensure mapping is cleared if invalid
                    index_to_heap_pos_[pos_index] = INVALID_POS;
                 }
                 const PositionType new_heap_pos = heap_data_.size();
                 heap_data_.emplace_back(std::move(value), index);
                 set_pos_unsafe(index, new_heap_pos);
                 sift_up(new_heap_pos);
            }
        }

        /**
         * @brief Inserts a new element. Assumes the index does not already exist.
         * If the index might exist, use insert_or_update or check contains() first.
         *
         * @param value The value of the element.
         * @param index The index of the element.
         */
        void insert(ValueType value, IndexType index)
            noexcept(noexcept(ensure_index_capacity(index)) &&
                     noexcept(heap_data_.emplace_back(std::move(value), index)) &&
                     noexcept(sift_up(0)))
        {
            // Assertion helps catch incorrect usage, but remove for release builds if needed.
            assert(!contains(index) && "Index already exists in PriorityQueue::insert. Use insert_or_update or remove first.");

            ensure_index_capacity(index);

            const PositionType new_pos = heap_data_.size();
            heap_data_.emplace_back(std::move(value), index);
            set_pos_unsafe(index, new_pos);
            sift_up(new_pos);
        }


        /**
         * @brief Decreases the key (value) of an existing element.
         * Does nothing if the index doesn't exist or if the new value is not smaller.
         *
         * @param new_value The new, smaller value.
         * @param index The index of the element to update.
         */
        void decrease_key(ValueType new_value, IndexType index)
            noexcept(noexcept(heap_comp_.value_comp(new_value, heap_data_[ROOT_POS].first))
                && noexcept(std::move(new_value)) && noexcept(sift_up(ROOT_POS)))
        {
            const PositionType pos_index = to_pos_index(index);

            // Check if index is potentially valid and mapped
            if (pos_index >= index_to_heap_pos_.size()) {
                return; // Index out of bounds for map
            }
            const PositionType pos = get_pos_unsafe(index); // Get potential position from map

            // Check if position is valid within heap and consistent
            if (!is_valid_and_consistent(pos, pos_index)) {
                 return; // Element not actually in the heap at the mapped position
            }

            // Check if the new value is actually smaller according to the comparator
            if (!heap_comp_.value_comp(new_value, heap_data_[pos].first)) {
                 // If new_value is not strictly smaller (could be equal or greater), do nothing.
                 // Original code had a check for equality, allowing update if equal. Keep that?
                 // Let's match original: if equal, update but don't sift.
                 if (heap_comp_.value_comp(heap_data_[pos].first, new_value)) {
                      // new_value is strictly greater, this is decrease_key, so ignore.
                     return;
                 }
                 // Here: new_value is equal or incomparable but not strictly smaller.
                 // Update if it's different (e.g. floating point non-equality)
                  if (!(new_value == heap_data_[pos].first)) { // Requires operator==
                     heap_data_[pos].first = std::move(new_value);
                  }
                 return; // Don't sift if not strictly smaller
            }

            // Update the value and sift up
            heap_data_[pos].first = std::move(new_value);
            sift_up(pos);
        }

        /**
         * @brief Deletes an element specified by its index.
         * Does nothing if the index does not exist in the heap.
         *
         * @param index The index of the element to delete.
         */
        void delete_element(IndexType index)
            noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_up(ROOT_POS)) && noexcept(sift_down(ROOT_POS)))
        {
            const PositionType pos_index = to_pos_index(index);

            // Check if index is potentially valid and mapped
            if (pos_index >= index_to_heap_pos_.size()) {
                return; // Index out of bounds for map
            }
            const PositionType pos = get_pos_unsafe(index); // Get potential position from map

            // Check if position is valid within heap and consistent
            if (!is_valid_and_consistent(pos, pos_index)) {
                return; // Element not actually in the heap at the mapped position
            }

            // Remove the element at the found position
            remove_at(pos);
        }
    };
} // namespace cluster_approx
