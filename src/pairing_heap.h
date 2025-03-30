#ifndef __PAIRING_HEAP_H__
#define __PAIRING_HEAP_H__

#include <vector>
#include <utility>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <stdexcept>

namespace cluster_approx {
    template <typename ValueType, typename PayloadType>
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
                    
                    return recycled_node;
                } else {
                    return new Node{nullptr, nullptr, nullptr, value, ValueType{}, std::move(payload), nullptr};
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
            
            explicit PairingHeap(std::vector<ItemHandle>* shared_buffer) : root(nullptr), buffer(shared_buffer), free_list_head(nullptr) {
                assert(buffer != nullptr && "Shared buffer cannot be null");
            }

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
                other.buffer = nullptr;
                other.free_list_head = nullptr;
            }

            PairingHeap& operator=(PairingHeap&& other) noexcept {
                if (this != &other) {
                    release_memory();

                    root = other.root;
                    buffer = other.buffer;
                    free_list_head = other.free_list_head;
                    other.root = nullptr;
                    other.buffer = nullptr;
                    other.free_list_head = nullptr;
                }
        
                return *this;
            }

            void release_memory() {
                if (root) {
                    std::vector<Node*> nodes_to_process_local;
                    std::vector<Node*>* nodes_to_process_ptr;
                    bool used_shared_buffer = false;

                    if (buffer && buffer->empty()) {
                        buffer->clear();
                        
                        nodes_to_process_ptr = buffer;
                        used_shared_buffer = true;
                    } else {
                        nodes_to_process_ptr = &nodes_to_process_local;
                
                        nodes_to_process_ptr->reserve(128);
                    }

                    std::vector<Node*>& nodes_to_process = *nodes_to_process_ptr;

                    nodes_to_process.push_back(root);
            
                    size_t current_idx = 0;
            
                    while(current_idx < nodes_to_process.size()) {
                        Node* node = nodes_to_process[current_idx++];
                        Node* child = node->child;
                
                        while(child) {
                            nodes_to_process.push_back(child);
                    
                            child = child->sibling;
                        }
                    }

                    for(Node* node_to_delete : nodes_to_process) {
                        delete node_to_delete;
                    }

                    if (used_shared_buffer) {
                        buffer->clear();
                    }
                }
                
                root = nullptr;

                clear_free_list();
            }

            inline bool is_empty() const noexcept {
                return root == nullptr;
            }

            inline bool get_min(ValueType* value, PayloadType* payload) const {
                if (root) {
                    assert(value != nullptr && payload != nullptr);
        
                    *value = root->value;
                    *payload = root->payload;
        
                    return true;
                }
        
                return false;
            }

            ItemHandle insert(ValueType value, PayloadType payload) {
                Node* new_node = allocate_node(value, std::move(payload));
                root = link(root, new_node);
        
                return new_node;
            }

            inline void add_to_heap(ValueType value) noexcept {
                if (root) {
                    root->value += value;
                    root->child_offset += value;
                }
            }

            void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
                assert(node != nullptr && "Cannot decrease key on a null handle");
                assert(!is_empty() && "Cannot decrease key on an empty heap");
                assert(to_value <= from_value && "Key decrease operation increased the key");

                ValueType stored_value_before_update = node->value;
                ValueType additional_offset = from_value - stored_value_before_update;
                node->child_offset += additional_offset;
                node->value = to_value;

                if (node == root) {
                    return;
                }

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
                root = link(root, node);
            }

            bool delete_min(ValueType* value, PayloadType* payload) {
                if (!root) {
                    return false;
                }
        
                assert(buffer != nullptr && "Cannot delete min without a valid buffer");
                assert(value != nullptr && payload != nullptr);

                Node* old_root = root;
                *value = old_root->value;
                *payload = std::move(old_root->payload);

                buffer->clear();

                Node* current_child = old_root->child;
        
                while (current_child) {
                    Node* next_sibling = current_child->sibling;
                    current_child->value += old_root->child_offset;
                    current_child->child_offset += old_root->child_offset;
                    current_child->left_up = nullptr;
                    current_child->sibling = nullptr;

                    buffer->push_back(current_child);
                    
                    current_child = next_sibling;
                }

                deallocate_node(old_root);
                
                root = nullptr;
                const size_t num_sub_heaps = buffer->size();
        
                if (num_sub_heaps == 0) {
                    return true;
                }

                if (num_sub_heaps == 1) {
                    root = (*buffer)[0];
                } else {
                    size_t merge_end = 0;

                    for (size_t i = 0; i + 1 < num_sub_heaps; i += 2) {
                        (*buffer)[merge_end++] = link((*buffer)[i], (*buffer)[i + 1]);
                    }
            
                    if (num_sub_heaps % 2 != 0) {
                        (*buffer)[merge_end++] = (*buffer)[num_sub_heaps - 1];
                    }

                    root = (*buffer)[merge_end - 1];
                    
                    for (size_t i = merge_end - 1; i > 0; --i) {
                        root = link((*buffer)[i - 1], root);
                    }
                }

                buffer->clear();
        
                return true;
            }

            static PairingHeap meld(PairingHeap* heap1, PairingHeap* heap2) {
                assert(heap1 != nullptr && heap2 != nullptr);
                assert(heap1->buffer != nullptr && heap2->buffer != nullptr);
                assert(heap1->buffer == heap2->buffer && "Heaps must share the same buffer for melding");

                PairingHeap result(heap1->buffer);

                result.root = link(heap1->root, heap2->root);

                if (heap1->free_list_head) {
                    Node* tail = heap1->free_list_head;
                    while (tail->next_free) {
                        tail = tail->next_free;
                    }
            
                    tail->next_free = heap2->free_list_head;
                    result.free_list_head = heap1->free_list_head;
                } else {
                    result.free_list_head = heap2->free_list_head;
                }

                heap1->root = nullptr;
                heap2->root = nullptr;
                heap1->free_list_head = nullptr;
                heap2->free_list_head = nullptr;

                return result;
            }

        private:
            Node* root;
            std::vector<ItemHandle>* buffer;
            Node* free_list_head;

            inline static Node* link(Node* node1, Node* node2) {
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

}

#endif
