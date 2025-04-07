#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <functional>
#include <type_traits>
#include <cstddef>
#include <algorithm>
#include <concepts>
#include <optional>

namespace cluster_approx {
    template<typename T>
    concept Indexable = std::is_integral_v<T> || std::is_enum_v<T>;

    namespace internal {
        template <typename T, bool = std::is_enum_v<T>>
        struct get_underlying_or_self_impl { using type = T; };

        template <typename T>
        struct get_underlying_or_self_impl<T, true> { using type = std::underlying_type_t<T>; };

        template <Indexable T>
        using get_underlying_or_self_t = typename get_underlying_or_self_impl<T>::type;

        template <typename V, Indexable I, typename C>
            requires std::strict_weak_order<C, V, V>
        struct HeapElementCompare {
            [[no_unique_address]] C value_comp;
            std::less<I> index_comp {};

            explicit HeapElementCompare(const C& comp = C()) noexcept(std::is_nothrow_copy_constructible_v<C>)
                : value_comp(comp) {}

            [[nodiscard]] constexpr bool operator()(const std::pair<V, I>& lhs, const std::pair<V, I>& rhs) const
                noexcept(noexcept(value_comp(lhs.first, rhs.first)) && noexcept(index_comp(lhs.second, rhs.second)))
            {
                if (value_comp(lhs.first, rhs.first)) return true;
                if (value_comp(rhs.first, lhs.first)) return false;

                return index_comp(lhs.second, rhs.second);
            }
        };
    }

    template <typename ValueType,
              Indexable IndexType,
              typename Compare = std::less<ValueType>>
        requires std::strict_weak_order<Compare, ValueType, ValueType>
              && std::totally_ordered<ValueType>
              && std::movable<ValueType>
              && std::movable<IndexType>
    class PriorityQueue {
        private:
            using PositionType = std::size_t;
            using UnderlyingIndexType = internal::get_underlying_or_self_t<IndexType>;
            using HeapElement = std::pair<ValueType, IndexType>;
            using Comp = internal::HeapElementCompare<ValueType, IndexType, Compare>;
            static constexpr PositionType INVALID_POS = std::numeric_limits<PositionType>::max();
            static constexpr PositionType ROOT_POS = 0;
            std::vector<HeapElement> heap_data_;
            std::vector<PositionType> index_to_heap_pos_;
            [[no_unique_address]] Comp heap_comp_;

            [[nodiscard]] static constexpr PositionType parent(PositionType i) noexcept { return (i - 1) / 2; }
            [[nodiscard]] static constexpr PositionType left(PositionType i) noexcept { return (i * 2) + 1; }
            [[nodiscard]] static constexpr PositionType right(PositionType i) noexcept { return (i * 2) + 2; }
            [[nodiscard]] static constexpr PositionType to_pos_index(IndexType index) noexcept {
                UnderlyingIndexType underlying_idx = static_cast<UnderlyingIndexType>(index);
                return static_cast<PositionType>(underlying_idx);
            }

            void ensure_index_capacity(IndexType index) {
                const PositionType pos_index = to_pos_index(index);
                if (pos_index >= index_to_heap_pos_.size()) {
                    PositionType current_size = index_to_heap_pos_.size();
                    PositionType desired_size = std::max(current_size + (current_size >> 1) + 1, pos_index + 1);
                    index_to_heap_pos_.resize(desired_size, INVALID_POS);
                }
            }

            constexpr void set_pos_unsafe(IndexType index, PositionType heap_pos) noexcept {
                index_to_heap_pos_[to_pos_index(index)] = heap_pos;
            }

            [[nodiscard]] constexpr PositionType get_pos_unsafe(IndexType index) const noexcept {
                return index_to_heap_pos_[to_pos_index(index)];
            }

            [[nodiscard]] constexpr bool is_valid_and_consistent(PositionType pos, PositionType pos_index) const noexcept {
                return pos != INVALID_POS
                    && pos < heap_data_.size()
                    && pos_index < index_to_heap_pos_.size()
                    && to_pos_index(heap_data_[pos].second) == pos_index;
            }

            [[nodiscard]] bool contains(IndexType index) const noexcept {
                const PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos_.size()) {
                    return false;
                }

                const PositionType heap_pos = index_to_heap_pos_[pos_index];

                return is_valid_and_consistent(heap_pos, pos_index);
            }

            void sift_up(PositionType i) noexcept(
                noexcept(std::move(heap_data_[i])) &&
                noexcept(heap_comp_(std::declval<HeapElement&>(), std::declval<HeapElement&>())) &&
                noexcept(set_pos_unsafe(std::declval<IndexType&>(), std::declval<PositionType&>()))
            )
            {
                if (i == ROOT_POS || i >= heap_data_.size()) return;

                HeapElement target = std::move(heap_data_[i]);
                PositionType current_pos = i;

                while (current_pos > ROOT_POS) {
                    PositionType p_pos = parent(current_pos);

                    if (!heap_comp_(target, heap_data_[p_pos])) {
                        break;
                    }

                    heap_data_[current_pos] = std::move(heap_data_[p_pos]);

                    set_pos_unsafe(heap_data_[current_pos].second, current_pos);

                    current_pos = p_pos;
                }

                heap_data_[current_pos] = std::move(target);

                set_pos_unsafe(heap_data_[current_pos].second, current_pos);
            }

            void sift_down(PositionType i) noexcept(
                noexcept(std::move(heap_data_[i])) &&
                noexcept(heap_comp_(std::declval<HeapElement&>(), std::declval<HeapElement&>())) &&
                noexcept(set_pos_unsafe(std::declval<IndexType&>(), std::declval<PositionType&>()))
            )
            {
                const PositionType heap_size = heap_data_.size();

                if (heap_size <= 1 || i >= heap_size) return;

                HeapElement target = std::move(heap_data_[i]);
                PositionType current_pos = i;

                while (true) {
                    PositionType child_pos = left(current_pos);

                    if (child_pos >= heap_size) {
                        break;
                    }

                    PositionType smaller_child_pos = find_smaller_child(current_pos, heap_size);

                    if (heap_comp_(heap_data_[smaller_child_pos], target)) {
                        heap_data_[current_pos] = std::move(heap_data_[smaller_child_pos]);

                        set_pos_unsafe(heap_data_[current_pos].second, current_pos);

                        current_pos = smaller_child_pos;
                    } else {
                        break;
                    }
                }

                heap_data_[current_pos] = std::move(target);

                set_pos_unsafe(heap_data_[current_pos].second, current_pos);
            }

            [[nodiscard]] PositionType find_smaller_child(PositionType i, PositionType heap_size) const noexcept {
                PositionType l_pos = left(i);
                PositionType r_pos = right(i);
                PositionType smaller_child_pos = l_pos;

                if (r_pos < heap_size && heap_comp_(heap_data_[r_pos], heap_data_[l_pos])) {
                    smaller_child_pos = r_pos;
                }

                return smaller_child_pos;
            }

            void remove_at(PositionType pos) noexcept (
                noexcept(std::move(heap_data_[pos])) &&
                noexcept(sift_up(pos)) &&
                noexcept(sift_down(pos))
                )
            {
                const PositionType pos_index = to_pos_index(heap_data_[pos].second);
                index_to_heap_pos_[pos_index] = INVALID_POS;
                const PositionType last_pos = heap_data_.size() - 1;

                if (pos == last_pos) {
                    heap_data_.pop_back();
                } else {
                    heap_data_[pos] = std::move(heap_data_[last_pos]);

                    heap_data_.pop_back();
                    set_pos_unsafe(heap_data_[pos].second, pos);

                    if (pos > ROOT_POS && heap_comp_(heap_data_[pos], heap_data_[parent(pos)])) {
                        sift_up(pos);
                    } else {
                        sift_down(pos);
                    }
                }
            }

        public:
            explicit PriorityQueue(const Compare& compare = Compare()) noexcept(noexcept(Comp(compare)))
                : heap_comp_(compare) {}

            PriorityQueue(PriorityQueue&&) noexcept = default;
            PriorityQueue& operator=(PriorityQueue&&) noexcept = default;
            PriorityQueue(const PriorityQueue&) = delete;
            PriorityQueue& operator=(const PriorityQueue&) = delete;
            ~PriorityQueue() = default;

            [[nodiscard]] bool is_empty() const noexcept {
                return heap_data_.empty();
            }

            [[nodiscard]] std::size_t size() const noexcept {
                return heap_data_.size();
            }

            [[nodiscard]] bool get_min(ValueType* value, IndexType* index) const
                noexcept(std::is_nothrow_copy_assignable_v<ValueType> && std::is_nothrow_copy_assignable_v<IndexType>)
            {
                if (is_empty()) {
                    if (value) *value = ValueType{};
                    if (index) *index = IndexType{};

                    return false;
                }
                if (!value || !index) {
                    return false;
                }

                *value = heap_data_[ROOT_POS].first;
                *index = heap_data_[ROOT_POS].second;

                return true;
            }

            [[nodiscard]] std::optional<HeapElement> get_min() const {
                if(is_empty()) {
                    return std::nullopt;
                }

                return heap_data_[ROOT_POS];
            }


            [[nodiscard]] bool delete_min(ValueType* value, IndexType* index)
                noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_down(ROOT_POS)))
            {
                if (is_empty()) {
                    if (value) *value = ValueType{};
                    if (index) *index = IndexType{};

                    return false;
                }
                if (!value || !index) {
                    return false;
                }

                *value = std::move(heap_data_[ROOT_POS].first);
                *index = heap_data_[ROOT_POS].second;

                remove_at(ROOT_POS);

                return true;
            }

            std::optional<HeapElement> delete_min()
                noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_down(ROOT_POS)))
            {
                if (is_empty()) {
                    return std::nullopt;
                }

                HeapElement min_element = std::move(heap_data_[ROOT_POS]);

                remove_at(ROOT_POS);

                return min_element;
            }


            void insert_or_update(ValueType value, IndexType index)
                noexcept(noexcept(ensure_index_capacity(index)) &&
                        noexcept(delete_element(index)) &&
                        noexcept(heap_data_.emplace_back(std::move(value), index)) &&
                        noexcept(sift_up(0)) && noexcept(sift_down(0)))
            {
                ensure_index_capacity(index);

                const PositionType pos_index = to_pos_index(index);
                const PositionType current_heap_pos = index_to_heap_pos_[pos_index];

                if (is_valid_and_consistent(current_heap_pos, pos_index)) {
                    ValueType& current_value = heap_data_[current_heap_pos].first;
                    bool value_changed = !(current_value == value);
                    bool decreased = value_changed && heap_comp_.value_comp(value, current_value);
                    bool increased = value_changed && heap_comp_.value_comp(current_value, value);

                    if (value_changed) {
                        heap_data_[current_heap_pos].first = std::move(value);

                        if (decreased) {
                            sift_up(current_heap_pos);
                        } else if (increased) {
                            sift_down(current_heap_pos);
                        }
                    }
                } else {
                    if (pos_index < index_to_heap_pos_.size()) {
                        index_to_heap_pos_[pos_index] = INVALID_POS;
                    }

                    const PositionType new_heap_pos = heap_data_.size();

                    heap_data_.emplace_back(std::move(value), index);
                    set_pos_unsafe(index, new_heap_pos);
                    sift_up(new_heap_pos);
                }
            }

            void insert(ValueType value, IndexType index)
                noexcept(noexcept(ensure_index_capacity(index)) &&
                        noexcept(heap_data_.emplace_back(std::move(value), index)) &&
                        noexcept(sift_up(0)))
            {
                ensure_index_capacity(index);

                const PositionType new_pos = heap_data_.size();

                heap_data_.emplace_back(std::move(value), index);
                set_pos_unsafe(index, new_pos);
                sift_up(new_pos);
            }

            void decrease_key(ValueType new_value, IndexType index)
                noexcept(noexcept(heap_comp_.value_comp(new_value, heap_data_[ROOT_POS].first))
                    && noexcept(std::move(new_value)) && noexcept(sift_up(ROOT_POS)))
            {
                const PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos_.size()) {
                    return;
                }

                const PositionType pos = get_pos_unsafe(index);

                if (!is_valid_and_consistent(pos, pos_index)) {
                    return;
                }

                if (!heap_comp_.value_comp(new_value, heap_data_[pos].first)) {
                    if (heap_comp_.value_comp(heap_data_[pos].first, new_value)) {
                        return;
                    }

                    if (!(new_value == heap_data_[pos].first)) {
                        heap_data_[pos].first = std::move(new_value);
                    }

                    return;
                }

                heap_data_[pos].first = std::move(new_value);

                sift_up(pos);
            }

            void delete_element(IndexType index)
                noexcept(noexcept(std::move(heap_data_[ROOT_POS])) && noexcept(sift_up(ROOT_POS)) && noexcept(sift_down(ROOT_POS)))
            {
                const PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos_.size()) {
                    return;
                }

                const PositionType pos = get_pos_unsafe(index);

                if (!is_valid_and_consistent(pos, pos_index)) {
                    return;
                }

                remove_at(pos);
            }
    };
}
