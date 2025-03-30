#pragma once

#include <vector>
#include <utility>
#include <optional>
#include <cstddef>
#include <cassert>
#include <concepts>

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
            Node* sibling;
            Node* child;
            Node* left_up;
            ValueType value;
            ValueType child_offset;
            PayloadType payload;
            Node* next_free;
        };

        Node* allocate_node(ValueType value, PayloadType payload) {
            if (free_list_head) {
                Node* recycled_node = free_list_head;
                free_list_head = recycled_node->next_free;
                recycled_node->sibling = nullptr;
                recycled_node->child = nullptr;
                recycled_node->left_up = nullptr;
                recycled_node->value = value;
                recycled_node->child_offset = ValueType{};
                recycled_node->payload = std::move(payload);
                recycled_node->next_free = nullptr;

                return recycled_node;
            } else {
                return new Node{
                    .sibling = nullptr,
                    .child = nullptr,
                    .left_up = nullptr,
                    .value = value,
                    .child_offset = ValueType{},
                    .payload = std::move(payload),
                    .next_free = nullptr
                };
            }
        }

        void deallocate_node(Node* node) {
            if (node) {
                node->next_free = free_list_head;
                free_list_head = node;
            }
        }

        void clear_free_list() {
            Node* current = free_list_head;
            
            while (current) {
                Node* next = current->next_free;
            
                delete current;
            
                current = next;
            }

            free_list_head = nullptr;
        }

    public:
        using ItemHandle = Node*;

        explicit PairingHeap(std::vector<ItemHandle>& shared_buffer)
            : root(nullptr), buffer(shared_buffer), free_list_head(nullptr) {}

        ~PairingHeap() {
            release_memory();
        }

        PairingHeap(const PairingHeap&) = delete;
        PairingHeap& operator=(const PairingHeap&) = delete;

        PairingHeap(PairingHeap&& other) noexcept
            : root(other.root),
              buffer(other.buffer),
              free_list_head(other.free_list_head)
        {
            other.root = nullptr;
            other.free_list_head = nullptr;
        }

        PairingHeap& operator=(PairingHeap&& other) noexcept {
            if (this != &other) {
                release_memory();

                root = other.root;
                buffer = other.buffer;
                free_list_head = other.free_list_head;
                other.root = nullptr;
                other.free_list_head = nullptr;
            }
            
            return *this;
        }

        void release_memory() {
            if (root) {
                buffer.clear();
                buffer.push_back(root);

                size_t current_idx = 0;
            
                while(current_idx < buffer.size()) {
                    Node* node = buffer[current_idx++];
                    Node* child = node->child;
            
                    while(child) {
                        buffer.push_back(child);
                        child = child->sibling;
                    }
                }

                for(Node* node_to_delete : buffer) {
                    delete node_to_delete;
                }

                buffer.clear();
            }
            
            root = nullptr;

            clear_free_list();
        }

        [[nodiscard]] bool is_empty() const noexcept {
            return root == nullptr;
        }

        [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
            if (root) {
                assert(value != nullptr && payload != nullptr && "Output pointers cannot be null");
            
                *value = root->value;
                *payload = root->payload;
            
                return true;
            }
            
            return false;
        }

        [[nodiscard]] std::optional<std::pair<ValueType, PayloadType>> get_min() const {
            if (root) {
                 return std::make_pair(root->value, root->payload);
            }
            
            return std::nullopt;
        }

        ItemHandle insert(ValueType value, PayloadType payload) {
            Node* new_node = allocate_node(value, std::move(payload));
            root = link(root, new_node);
            
            return new_node;
        }

        void add_to_heap(ValueType value) noexcept(noexcept(root->value += value) && noexcept(root->child_offset += value)) {
            if (root) {
                static_assert(noexcept(std::declval<ValueType&>() += std::declval<const ValueType&>()), "ValueType::operator+= must be noexcept for PairingHeap::add_to_heap to be noexcept");
                
                root->value += value;
                root->child_offset += value;
            }
        }

        void decrease_key(ItemHandle node, ValueType to_value) {
            assert(node != nullptr && "Cannot decrease key on a null handle");
            assert(!is_empty() && "Cannot decrease key on an empty heap");

            if (node == root) {
                assert(to_value <= node->value && "Cannot increase root key via decrease_key");
                node->value = to_value;
                
                return;
            }

            assert(node->left_up != nullptr && "Non-root node must have a left_up pointer");

            Node* parent_or_left_sibling = node->left_up;
            Node* right_sibling = node->sibling;

            assert(to_value <= node->value && "Decrease key 'to_value' must be <= node's internal value for standard pairing heap decrease-key.");
            
            node->value = to_value;

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
            root = link(root, node);
        }


        [[nodiscard]] bool delete_min(ValueType* value, PayloadType* payload) {
            if (!root) {
                return false;
            }
            
            assert(value != nullptr && payload != nullptr && "Output pointers cannot be null");

            Node* old_root = root;
            *value = old_root->value;
            *payload = std::move(old_root->payload);

            buffer.clear();

            Node* current_child = old_root->child;

            while (current_child) {
                Node* next_sibling = current_child->sibling;
                current_child->value += old_root->child_offset;
                current_child->child_offset += old_root->child_offset;
                current_child->left_up = nullptr;
                current_child->sibling = nullptr;

                buffer.push_back(current_child);
            
                current_child = next_sibling;
            }

            deallocate_node(old_root);

            root = nullptr;
            const size_t num_sub_heaps = buffer.size();
            
            if (num_sub_heaps == 0) {
                return true;
            }

            if (num_sub_heaps == 1) {
                root = buffer[0];
            } else {
                size_t merge_end = 0;
                
                for (size_t i = 0; i + 1 < num_sub_heaps; i += 2) {
                    buffer[merge_end++] = link(buffer[i], buffer[i + 1]);
                }

                if (num_sub_heaps % 2 != 0) {
                    buffer[merge_end++] = buffer[num_sub_heaps - 1];
                }

                if (merge_end > 0) {
                    root = buffer[merge_end - 1];
                    for (size_t i = merge_end - 1; i > 0; --i) {
                        root = link(buffer[i - 1], root);
                    }
                }

            }

            buffer.clear();

            return true;
        }

        [[nodiscard]] static PairingHeap meld(PairingHeap&& heap1, PairingHeap&& heap2) {
            assert(&heap1.buffer == &heap2.buffer && "Heaps must share the same buffer object for melding");

            PairingHeap result(heap1.buffer);
            result.root = link(heap1.root, heap2.root);

            if (heap1.free_list_head) {
                Node* tail = heap1.free_list_head;
                
                while (tail->next_free) {
                    tail = tail->next_free;
                }
                
                tail->next_free = heap2.free_list_head;
                result.free_list_head = heap1.free_list_head;
            } else {
                result.free_list_head = heap2.free_list_head;
            }

            heap1.root = nullptr;
            heap1.free_list_head = nullptr;
            heap2.root = nullptr;
            heap2.free_list_head = nullptr;

            return result;
        }


    private:
        Node* root;
        std::vector<ItemHandle>& buffer;
        Node* free_list_head;

        static Node* link(Node* node1, Node* node2) {
            if (!node2) return node1;
            if (!node1) return node2;

            Node* smaller_node;
            Node* larger_node;

            if (node1->value <= node2->value) {
                smaller_node = node1;
                larger_node = node2;
            } else {
                smaller_node = node2;
                larger_node = node1;
            }

            larger_node->value -= smaller_node->child_offset;
            larger_node->child_offset -= smaller_node->child_offset;
            larger_node->left_up = smaller_node;
            larger_node->sibling = smaller_node->child;

            if (smaller_node->child) {
                smaller_node->child->left_up = larger_node;
            }

            smaller_node->child = larger_node;

            return smaller_node;
        }
    };

} // namespace cluster_approx
