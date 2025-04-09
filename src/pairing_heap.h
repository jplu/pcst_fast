#pragma once

#include <vector>
#include <utility>
#include <optional>
#include <cstddef>
#include <cassert>
#include <concepts>
#include <stdexcept>
#include <algorithm>

namespace cluster_approx {
    template <typename T>
    concept PairingHeapValue = std::totally_ordered<T> && requires(T a, T b) {
        { a += b } -> std::same_as<T&>;
        { a -= b } -> std::same_as<T&>;
        T{};
    };

    template <PairingHeapValue ValueType, typename PayloadType>
    class PairingHeap {
    private:
        struct Node {
            Node* sibling = nullptr;
            Node* child = nullptr;
            Node* left_up = nullptr;
            ValueType value = {};
            ValueType child_offset = {};
            PayloadType payload = {};
            Node* next_free = nullptr;
        };

        Node* allocate_node(ValueType value, PayloadType payload) {
            if (free_list_head_) {
                return reuse_node(value, std::move(payload));
            } else {
                return new Node{
                    .value = value,
                    .payload = std::move(payload)
                };
            }
        }

        Node* reuse_node(ValueType value, PayloadType payload) {
            Node* recycled_node = free_list_head_;
            free_list_head_ = recycled_node->next_free;
            recycled_node->sibling = nullptr;
            recycled_node->child = nullptr;
            recycled_node->left_up = nullptr;
            recycled_node->value = value;
            recycled_node->child_offset = ValueType{};
            recycled_node->payload = std::move(payload);
            recycled_node->next_free = nullptr;
            return recycled_node;
        }

        void deallocate_node(Node* node) noexcept {
            if (node) {
                node->next_free = free_list_head_;
                free_list_head_ = node;
            }
        }

        void clear_free_list() noexcept {
            Node* current = free_list_head_;
            while (current) {
                Node* next = current->next_free;
                delete current;
                current = next;
            }
            free_list_head_ = nullptr;
        }

        void destroy_nodes(Node* node) noexcept {
            if (!node) return;
            std::vector<Node*> nodes_to_delete;
            nodes_to_delete.push_back(node);
            while(!nodes_to_delete.empty()) {
                Node* current = nodes_to_delete.back();
                nodes_to_delete.pop_back();
                Node* child = current->child;
                while (child) {
                    nodes_to_delete.push_back(child);
                    child = child->sibling;
                }
                delete current;
            }
        }

        static Node* link(Node* node1, Node* node2) noexcept {
             if (node1 == NULL) {
                 return node2;
             }
             if (node2 == NULL) {
                 return node1;
             }
             Node* smaller_node = node2;
             Node* larger_node = node1;
             // Compare presumably absolute values (as prepared by delete_min)
             if (node1->value < node2->value) {
                 smaller_node = node1;
                 larger_node = node2;
             }

             // *** Replicate EXACT old link structure update and adjustments ***
             larger_node->sibling = smaller_node->child;
             if (larger_node->sibling != NULL) {
                 larger_node->sibling->left_up = larger_node;
             }
             larger_node->left_up = smaller_node;
             smaller_node->child = larger_node;
             larger_node->value -= smaller_node->child_offset; // Old logic adjustment 1
             larger_node->child_offset -= smaller_node->child_offset; // Old logic adjustment 2

             return smaller_node;
         }

        Node* merge_sub_heaps(std::vector<Node*>& sub_heaps) noexcept {
            if (sub_heaps.empty()) {
                return nullptr;
            }
            if (sub_heaps.size() == 1) {
                return sub_heaps[0];
            }

            size_t merge_end = 0;
            for (size_t i = 0; i + 1 < sub_heaps.size(); i += 2) {
                sub_heaps[merge_end++] = link(sub_heaps[i], sub_heaps[i + 1]);
            }
            if (sub_heaps.size() % 2 != 0) {
                sub_heaps[merge_end++] = sub_heaps.back();
            }

            Node* merged_root = sub_heaps[merge_end - 1];
            for (size_t i = merge_end - 1; i > 0; --i) {
                merged_root = link(sub_heaps[i - 1], merged_root);
            }

            return merged_root;
        }

        void detach_node(Node* node) noexcept {
            assert(node != nullptr && node != root_ && "Cannot detach root or null node this way");
            assert(node->left_up != nullptr && "Non-root node must have a left_up pointer");

            Node* parent_or_left_sibling = node->left_up;
            Node* right_sibling = node->sibling;

            if (parent_or_left_sibling->child == node) {
                parent_or_left_sibling->child = right_sibling;
            } else {
                assert(parent_or_left_sibling->sibling == node && "Node's left_up inconsistent with sibling chain");
                parent_or_left_sibling->sibling = right_sibling;
            }

            if (right_sibling) {
                right_sibling->left_up = parent_or_left_sibling;
            }

            node->left_up = nullptr;
            node->sibling = nullptr;
        }

    public:
        using ItemHandle = Node*;
        static constexpr ItemHandle kInvalidHandle = nullptr;

        PairingHeap() = default;
        ~PairingHeap() noexcept {
            release_memory();
        }
        PairingHeap(const PairingHeap&) = delete;
        PairingHeap& operator=(const PairingHeap&) = delete;
        PairingHeap(PairingHeap&& other) noexcept
            : root_(other.root_),
              free_list_head_(other.free_list_head_)
        {
            other.root_ = nullptr;
            other.free_list_head_ = nullptr;
        }
        PairingHeap& operator=(PairingHeap&& other) noexcept {
            if (this != &other) {
                release_memory();
                root_ = other.root_;
                free_list_head_ = other.free_list_head_;
                other.root_ = nullptr;
                other.free_list_head_ = nullptr;
            }
            return *this;
        }

        void release_memory() noexcept {
            destroy_nodes(root_);
            root_ = nullptr;
            clear_free_list();
        }

        [[nodiscard]] bool is_empty() const noexcept {
            return root_ == nullptr;
        }

        [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
             if (!root_) return false;
             assert(value != nullptr && payload != nullptr && "Output pointers cannot be null");
             *value = root_->value;
             *payload = root_->payload;
             return true;
        }

        [[nodiscard]] std::optional<std::pair<ValueType, PayloadType>> get_min() const {
            if (root_) {
                 return std::make_pair(root_->value, root_->payload);
            }
            return std::nullopt;
        }

        ItemHandle insert(ValueType value, PayloadType payload) {
            Node* new_node = allocate_node(value, std::move(payload));
            root_ = link(root_, new_node);
            return new_node;
        }

        void add_to_heap(const ValueType& value) noexcept(noexcept(root_->value += value) && noexcept(root_->child_offset += value)) {
            static_assert(noexcept(std::declval<ValueType&>() += std::declval<const ValueType&>()),
                          "ValueType::operator+= must be noexcept for PairingHeap::add_to_heap to be fully noexcept");

            if (root_) {
                root_->value += value;
                root_->child_offset += value;
            }
        }

        void decrease_key(ItemHandle node, const ValueType& new_value) {
            if (!node) {
                throw std::invalid_argument("PairingHeap::decrease_key: Cannot operate on a null handle");
            }

            // Handle root node separately (its value is absolute)
            if (node == root_) {
                if (new_value > root_->value) {
                     // Allow equality, prevent increase
                     throw std::invalid_argument("PairingHeap::decrease_key: New value must be <= current value for root.");
                }
                root_->value = new_value;
                return;
            }

            // For non-root nodes:
            // We skip validation new_value <= current_absolute_value because it's hard to calculate.
            // Trust the caller provided a value representing a decrease.

            // Store the target absolute value directly, temporarily.
            node->value = new_value;

            // Detach the node from its current position and relink it to the root
            detach_node(node);
            root_ = link(root_, node);
            // Link is expected to make node->value relative again if node becomes a child.
        }

        [[nodiscard]] bool delete_min(ValueType* value, PayloadType* payload) {
            if (!root_) {
                return false;
            }
            assert(value != nullptr && payload != nullptr && "Output pointers cannot be null");

            Node* old_root = root_;
            *value = old_root->value; // Root value is absolute
            *payload = std::move(old_root->payload);
            std::vector<Node*> sub_heaps;
            Node* current_child = old_root->child;

            while (current_child) {
                Node* next_sibling = current_child->sibling;
                // Restore: Make child's value absolute relative to root for consistent linking.
                current_child->value += old_root->value; // Use root's absolute value
                // Restore: Propagate the offset change.
                current_child->child_offset += old_root->child_offset; // RESTORED
                current_child->left_up = nullptr;
                current_child->sibling = nullptr;
                sub_heaps.push_back(current_child);
                current_child = next_sibling;
            }

            deallocate_node(old_root);
            // merge_sub_heaps uses link. Link expects comparable (absolute) values
            // and makes the child relative using the parent's offset info.
            root_ = merge_sub_heaps(sub_heaps);
            return true;
        }

        [[nodiscard]] static PairingHeap meld(PairingHeap&& heap1, PairingHeap&& heap2) noexcept {
            PairingHeap result;
            result.root_ = link(heap1.root_, heap2.root_);

            // Combine free lists efficiently
            if (heap1.free_list_head_) {
                if (heap2.free_list_head_) {
                    // Find tail of heap1's free list
                    Node* tail = heap1.free_list_head_;
                    while (tail->next_free) {
                        tail = tail->next_free;
                    }
                    // Append heap2's list
                    tail->next_free = heap2.free_list_head_;
                }
                result.free_list_head_ = heap1.free_list_head_;
            } else {
                result.free_list_head_ = heap2.free_list_head_;
            }

            heap1.root_ = nullptr;
            heap1.free_list_head_ = nullptr;
            heap2.root_ = nullptr;
            heap2.free_list_head_ = nullptr;

            return result;
        }

    private:
        Node* root_ = nullptr;
        Node* free_list_head_ = nullptr;
    };
}