#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <functional>
#include <type_traits>
#include <cstddef>
#include <algorithm>
#include <concepts>


namespace cluster_approx {
    template<typename T>
    concept Indexable = std::is_integral_v<T> || std::is_enum_v<T>;

    namespace internal {
        template <typename T, bool = std::is_enum_v<T>>
        struct get_underlying_or_self_impl {
            using type = T;
        };

        template <typename T>
        struct get_underlying_or_self_impl<T, true> {
            using type = std::underlying_type_t<T>;
        };

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
            std::vector<HeapElement> heap_data_;
            std::vector<PositionType> index_to_heap_pos_;
            [[no_unique_address]] Comp heap_comp_;

            [[nodiscard]] static constexpr PositionType parent(PositionType i) noexcept { return (i - 1) / 2; }
            [[nodiscard]] static constexpr PositionType left(PositionType i) noexcept { return (i * 2) + 1; }
            [[nodiscard]] static constexpr PositionType right(PositionType i) noexcept { return (i * 2) + 2; }
            [[nodiscard]] static constexpr PositionType to_pos_index(IndexType index) noexcept {
                static_assert(std::is_unsigned_v<UnderlyingIndexType> || std::is_signed_v<UnderlyingIndexType>,
                            "Underlying index type must be integral");
                if constexpr (std::is_signed_v<UnderlyingIndexType>) {}

                return static_cast<PositionType>(static_cast<UnderlyingIndexType>(index));
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

            void sift_up(PositionType i) noexcept(
                noexcept(std::move(heap_data_[i])) &&
                noexcept(heap_comp_(std::declval<HeapElement&>(), std::declval<HeapElement&>())) &&
                noexcept(set_pos_unsafe(std::declval<IndexType&>(), std::declval<PositionType&>()))
            )
            {
                if (i == 0 || i >= heap_data_.size()) return;

                HeapElement target = std::move(heap_data_[i]);
                PositionType current_pos = i;

                while (current_pos > 0) {
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
                PositionType child_pos = left(current_pos);

                while (child_pos < heap_size) {
                    PositionType r_pos = right(current_pos);
                    PositionType smaller_child_pos = child_pos;

                    if (r_pos < heap_size && heap_comp_(heap_data_[r_pos], heap_data_[child_pos])) {
                        smaller_child_pos = r_pos;
                    }

                    if (heap_comp_(heap_data_[smaller_child_pos], target)) {
                        heap_data_[current_pos] = std::move(heap_data_[smaller_child_pos]);
                
                        set_pos_unsafe(heap_data_[current_pos].second, current_pos);
                
                        current_pos = smaller_child_pos;
                        child_pos = left(current_pos);
                    } else {
                        break;
                    }
                }
                
                heap_data_[current_pos] = std::move(target);
                
                set_pos_unsafe(heap_data_[current_pos].second, current_pos);
            }

        public:
            explicit PriorityQueue(const Compare& compare = Compare()) noexcept(noexcept(Comp(compare)))
                : heap_comp_(compare) {}

            PriorityQueue(const PriorityQueue&) = delete;
            PriorityQueue& operator=(const PriorityQueue&) = delete;
            PriorityQueue(PriorityQueue&&) noexcept = default;
            PriorityQueue& operator=(PriorityQueue&&) noexcept = default;
            ~PriorityQueue() = default;

            [[nodiscard]] bool is_empty() const noexcept {
                return heap_data_.empty();
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

                *value = heap_data_[0].first;
                *index = heap_data_[0].second;

                return true;
            }

            [[nodiscard]] bool delete_min(ValueType* value, IndexType* index)
                noexcept(noexcept(std::move(heap_data_[0])) && noexcept(sift_down(0)))
            {
                if (is_empty()) {
                    if (value) *value = ValueType{};
                    if (index) *index = IndexType{};
             
                    return false;
                }
                if (!value || !index) {
                    return false;
                }

                *value = std::move(heap_data_[0].first);
                *index = heap_data_[0].second;
                const PositionType pos_index = to_pos_index(*index);
                
                if (pos_index < index_to_heap_pos_.size()) {
                    index_to_heap_pos_[pos_index] = INVALID_POS;
                }

                const PositionType last_pos = heap_data_.size() - 1;
                
                if (last_pos > 0) {
                    heap_data_[0] = std::move(heap_data_[last_pos]);
                
                    set_pos_unsafe(heap_data_[0].second, 0);
                    heap_data_.pop_back();
                    sift_down(0);
                } else {
                    heap_data_.pop_back();
                }

                return true;
            }

            void insert(ValueType value, IndexType index)
                noexcept(noexcept(ensure_index_capacity(index)) && noexcept(delete_element(index))
                    && noexcept(heap_data_.emplace_back(std::move(value), index)) && noexcept(sift_up(0)))
            {
                ensure_index_capacity(index);
                
                const PositionType pos_index = to_pos_index(index);
                const PositionType existing_pos = (pos_index < index_to_heap_pos_.size()) ? index_to_heap_pos_[pos_index] : INVALID_POS;

                if (is_valid_and_consistent(existing_pos, pos_index)) {
                    delete_element(index);
                } else if (pos_index < index_to_heap_pos_.size()) {
                    index_to_heap_pos_[pos_index] = INVALID_POS;
                }

                const PositionType new_pos = heap_data_.size();
                
                heap_data_.emplace_back(std::move(value), index);
                set_pos_unsafe(index, new_pos);
                sift_up(new_pos);
            }

            void decrease_key(ValueType new_value, IndexType index)
                noexcept(noexcept(heap_comp_.value_comp(new_value, heap_data_[0].first))
                    && noexcept(std::move(new_value)) && noexcept(sift_up(0)))
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
                
                    heap_data_[pos].first = std::move(new_value);
                
                    return;
                }

                heap_data_[pos].first = std::move(new_value);
                
                sift_up(pos);
            }

            void delete_element(IndexType index)
                noexcept(noexcept(std::move(heap_data_[0])) && noexcept(sift_up(0)) && noexcept(sift_down(0)))
            {
                const PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos_.size()) {
                    return;
                }

                const PositionType pos = get_pos_unsafe(index);

                if (!is_valid_and_consistent(pos, pos_index)) {
                    return;
                }

                index_to_heap_pos_[pos_index] = INVALID_POS;
                const PositionType last_pos = heap_data_.size() - 1;

                if (pos == last_pos) {
                    heap_data_.pop_back();
                } else {
                    heap_data_[pos] = std::move(heap_data_[last_pos]);
                    
                    heap_data_.pop_back();

                    set_pos_unsafe(heap_data_[pos].second, pos);

                    if (pos > 0 && heap_comp_(heap_data_[pos], heap_data_[parent(pos)])) {
                        sift_up(pos);
                    } else {
                        sift_down(pos);
                    }
                }
            }
    };
}
