#ifndef __PRIORITY_QUEUE_H__
#define __PRIORITY_QUEUE_H__

#include <vector>
#include <utility>
#include <limits>
#include <functional>
#include <type_traits>
#include <cstddef>
#include <stdexcept>
#include <algorithm>


namespace cluster_approx {
    namespace internal {
        template <typename T, bool is_enum = std::is_enum<T>::value>
        struct get_underlying_or_self {
            using type = T;
        };

        template <typename T>
        struct get_underlying_or_self<T, true> {
            using type = typename std::underlying_type<T>::type;
        };

        template <typename V, typename I, typename C = std::less<V>>
        struct HeapElementCompare {
            C value_comp;
            std::less<I> index_comp;

            explicit HeapElementCompare(const C& comp = C()) noexcept(std::is_nothrow_default_constructible<C>::value && std::is_nothrow_copy_constructible<C>::value && std::is_nothrow_default_constructible<std::less<I>>::value)
            : value_comp(comp), index_comp() {}

            inline bool operator()(const std::pair<V, I>& lhs, const std::pair<V, I>& rhs) const {
                if (value_comp(lhs.first, rhs.first)) return true;
                if (value_comp(rhs.first, lhs.first)) return false;
                return index_comp(lhs.second, rhs.second);
            }
        };

    }

    template <typename ValueType,
            typename IndexType,
            typename Compare = std::less<ValueType>>
    class PriorityQueue {
        private:
            using HeapElement = std::pair<ValueType, IndexType>;
            using PositionType = std::size_t;
            static constexpr PositionType INVALID_POS = std::numeric_limits<PositionType>::max();
            using UnderlyingIndexType = typename internal::get_underlying_or_self<IndexType>::type;
            std::vector<HeapElement> heap_data;
            std::vector<PositionType> index_to_heap_pos;
            internal::HeapElementCompare<ValueType, IndexType, Compare> heap_comp;

            inline PositionType parent(PositionType i) const noexcept { return (i - 1) >> 1; }
            inline PositionType left(PositionType i) const noexcept { return (i << 1) + 1; }
            inline PositionType right(PositionType i) const noexcept { return (i << 1) + 2; }
            inline PositionType to_pos_index(IndexType index) const noexcept {
                return static_cast<PositionType>(static_cast<UnderlyingIndexType>(index));
            }

            void ensure_index_capacity(IndexType index) {
                PositionType pos_index = to_pos_index(index);
                if (pos_index >= index_to_heap_pos.size()) {
                    PositionType current_size = index_to_heap_pos.size();
                    PositionType desired_size = std::max(static_cast<PositionType>(current_size + (current_size >> 1)), pos_index + 1);

                    index_to_heap_pos.resize(desired_size, INVALID_POS);
                }
            }

            inline void set_pos_unsafe(IndexType index, PositionType heap_pos) noexcept {
                index_to_heap_pos[to_pos_index(index)] = heap_pos;
            }

            void sift_up(PositionType i) {
                HeapElement target = std::move(heap_data[i]);
                PositionType current_pos = i;

                while (current_pos > 0) {
                    PositionType p = parent(current_pos);
                    
                    if (!heap_comp(target, heap_data[p])) {
                        break;
                    }
                    
                    heap_data[current_pos] = std::move(heap_data[p]);
                    
                    set_pos_unsafe(heap_data[current_pos].second, current_pos);
                    
                    current_pos = p;
                }
                
                heap_data[current_pos] = std::move(target);
                
                set_pos_unsafe(heap_data[current_pos].second, current_pos);
            }


            void sift_down(PositionType i) {
                PositionType heap_size = heap_data.size();
                
                if (heap_size <= 1 || i >= heap_size) return;

                HeapElement target = std::move(heap_data[i]);
                PositionType current_pos = i;
                PositionType l = left(current_pos);

                while (l < heap_size) {
                    PositionType r = right(current_pos);
                    PositionType smallest_child_pos = l;

                    if (r < heap_size && heap_comp(heap_data[r], heap_data[l])) {
                        smallest_child_pos = r;
                    }

                    if (heap_comp(heap_data[smallest_child_pos], target)) {
                        heap_data[current_pos] = std::move(heap_data[smallest_child_pos]);
                        
                        set_pos_unsafe(heap_data[current_pos].second, current_pos);
                        
                        current_pos = smallest_child_pos; // Move hole down
                        l = left(current_pos);
                    } else {
                        break;
                    }
                }

                heap_data[current_pos] = std::move(target);
                
                set_pos_unsafe(heap_data[current_pos].second, current_pos);
            }


            inline PositionType get_pos_unsafe(IndexType index) const noexcept {
                return index_to_heap_pos[to_pos_index(index)];
            }

            inline bool is_valid_and_consistent(PositionType pos, PositionType pos_index) const noexcept {
                return pos != INVALID_POS && pos < heap_data.size() && to_pos_index(heap_data[pos].second) == pos_index;
            }


        public:
            explicit PriorityQueue(const Compare& compare = Compare())
            : heap_comp(compare) {}

            PriorityQueue(const PriorityQueue&) = delete;
            PriorityQueue& operator=(const PriorityQueue&) = delete;

            bool is_empty() const noexcept {
                return heap_data.empty();
            }

            bool get_min(ValueType* value, IndexType* index) const {
                if (is_empty() || !value || !index) {
                    return false;
                }

                *value = heap_data[0].first;
                *index = heap_data[0].second;

                return true;
            }

            bool delete_min(ValueType* value, IndexType* index) {
                if (is_empty() || !value || !index) {
                    return false;
                }

                *value = heap_data[0].first;
                *index = heap_data[0].second;
                PositionType pos_index = to_pos_index(*index);
                
                if (pos_index < index_to_heap_pos.size()) {
                    index_to_heap_pos[pos_index] = INVALID_POS;
                }

                PositionType last_pos = heap_data.size() - 1;
                
                if (last_pos > 0) {
                    heap_data[0] = std::move(heap_data[last_pos]);
                    PositionType moved_elem_pos_index = to_pos_index(heap_data[0].second);
                
                    if (moved_elem_pos_index < index_to_heap_pos.size()) {
                        set_pos_unsafe(heap_data[0].second, 0);
                    }

                    heap_data.pop_back();
                    sift_down(0);
                } else {
                    heap_data.pop_back();
                }
                
                return true;
            }

            void insert(ValueType value, IndexType index) {
                ensure_index_capacity(index);
                
                PositionType pos_index = to_pos_index(index);
                PositionType existing_pos = index_to_heap_pos[pos_index];
                
                if (is_valid_and_consistent(existing_pos, pos_index)) {
                    delete_element(index);
                } else if (existing_pos != INVALID_POS) {
                    index_to_heap_pos[pos_index] = INVALID_POS;
                }
                
                PositionType new_pos = heap_data.size();
                
                heap_data.emplace_back(std::move(value), index);
                set_pos_unsafe(index, new_pos);
                sift_up(new_pos);
            }

            void decrease_key(ValueType new_value, IndexType index) {
                PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos.size()) {
                    return;
                }

                PositionType pos = get_pos_unsafe(index);

                if (!is_valid_and_consistent(pos, pos_index)) {
                    return;
                }

                if (!heap_comp.value_comp(new_value, heap_data[pos].first)) {
                    return;
                }

                heap_data[pos].first = std::move(new_value);
                
                sift_up(pos);
            }

            void delete_element(IndexType index) {
                PositionType pos_index = to_pos_index(index);

                if (pos_index >= index_to_heap_pos.size()) {
                    return;
                }

                PositionType pos = get_pos_unsafe(index);

                if (!is_valid_and_consistent(pos, pos_index)) {
                    return;
                }

                PositionType last_pos = heap_data.size() - 1;
                index_to_heap_pos[pos_index] = INVALID_POS;

                if (pos == last_pos) {
                    heap_data.pop_back();
                } else {
                    heap_data[pos] = std::move(heap_data[last_pos]);
                    PositionType moved_elem_pos_index = to_pos_index(heap_data[pos].second);

                    if (moved_elem_pos_index < index_to_heap_pos.size()) {
                        set_pos_unsafe(heap_data[pos].second, pos);
                    }

                    heap_data.pop_back();

                    if (pos > 0 && heap_comp(heap_data[pos], heap_data[parent(pos)])) {
                        sift_up(pos);
                    } else {
                        sift_down(pos);
                    }
                }
            }
    };

    template <typename ValueType, typename IndexType, typename Compare>
    constexpr typename PriorityQueue<ValueType, IndexType, Compare>::PositionType
        PriorityQueue<ValueType, IndexType, Compare>::INVALID_POS;
}

#endif
