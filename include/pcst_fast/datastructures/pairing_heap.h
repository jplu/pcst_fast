#pragma once

#include <vector>
#include <cstdlib>
#include <utility>
#include <limits>
#include <cassert>

namespace cluster_approx {

/**
 * @brief Implements a Pairing Heap data structure.
 *
 * Supports insert, delete-min, decrease-key, and meld operations.
 * It uses lazy updates for additions to the entire heap and during linking.
 * Note: This implementation relies on an external shared buffer for temporary
 * storage during delete_min and release_memory operations to potentially
 * reduce allocations. The user must ensure the provided buffer pointer is valid
 * during the heap's lifetime.
 *
 * @tparam ValueType The type of the values (keys) stored in the heap. Must support comparison and arithmetic operations (+, -).
 * @tparam PayloadType The type of the payload associated with each value.
 */
template <typename ValueType, typename PayloadType>
class PairingHeap {
  private:
    struct Node {
        Node* sibling = nullptr;
        Node* child = nullptr;
        Node* left_up = nullptr;
        ValueType value = ValueType{};
        ValueType child_offset = ValueType{};
        PayloadType payload = PayloadType{};
    };

  public:
    using ItemHandle = Node*;

    /**
     * @brief Constructs a PairingHeap.
     * @param shared_buffer A pointer to a vector managed externally, used for temporary storage. Must not be null.
     */
    explicit PairingHeap(std::vector<ItemHandle>* shared_buffer) : root_(nullptr), buffer_(shared_buffer) {
        assert(buffer_ != nullptr && "Shared buffer cannot be null.");
    }

    PairingHeap(const PairingHeap&) = delete;
    PairingHeap& operator=(const PairingHeap&) = delete;

    PairingHeap(PairingHeap&& other) noexcept
        : root_(other.root_), buffer_(other.buffer_) {
        other.root_ = nullptr;

    }

    PairingHeap& operator=(PairingHeap&& other) noexcept {
        if (this != &other) {
            release_memory();
            root_ = other.root_;
            buffer_ = other.buffer_;
            other.root_ = nullptr;
        }
        return *this;
    }

    /**
     * @brief Releases all memory allocated for the heap nodes.
     * Requires the shared buffer passed during construction.
     */
    void release_memory() {
        if (root_ == nullptr) {
            return;
        }
        assert(buffer_ != nullptr && "Shared buffer required for release_memory.");
        buffer_->clear();
        buffer_->push_back(root_);
        size_t current_idx = 0;
        while (current_idx < buffer_->size()) {
            Node* current_node = (*buffer_)[current_idx];
            assert(current_node != nullptr);
            if (current_node->child != nullptr) {
                buffer_->push_back(current_node->child);
            }
            if (current_node->sibling != nullptr) {
                buffer_->push_back(current_node->sibling);
            }
            current_idx++;
        }

        for (Node* node_to_delete : *buffer_) {
            delete node_to_delete;
        }
        buffer_->clear();
        root_ = nullptr;
    }

    /**
     * @brief Checks if the heap is empty.
     * @return True if the heap contains no elements, false otherwise.
     */
    [[nodiscard]] bool is_empty() const noexcept {
        return root_ == nullptr;
    }

    /**
     * @brief Gets the minimum value and payload without removing the element.
     * @param[out] value Pointer to store the minimum value.
     * @param[out] payload Pointer to store the payload of the minimum element.
     * @return True if the heap is not empty and values were retrieved, false otherwise.
     */
    [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
        assert(value != nullptr && payload != nullptr);
        if (root_ != nullptr) {
            *value = root_->value;
            *payload = root_->payload;
            return true;
        } else {
            return false;
        }
    }

    /**
     * @brief Inserts a new element into the heap.
     * @param value The value (key) of the element.
     * @param payload The payload associated with the value.
     * @return An ItemHandle representing the newly inserted node.
     */
    [[nodiscard]] ItemHandle insert(ValueType value, PayloadType payload) {
        Node* new_node = new Node();

        new_node->value = value;
        new_node->payload = payload;
        root_ = link(root_, new_node);
        assert(root_ != nullptr);
        return new_node;
    }

    /**
     * @brief Adds a value to all elements currently in the heap (lazy update).
     * This is implemented by adding the value to the root's value and offset.
     * @param value The value to add.
     */
    void add_to_heap(ValueType value) {
        if (root_ != nullptr) {
            root_->value += value;
            root_->child_offset += value;
        }
    }

    /**
     * @brief Decreases the key (value) of a specific node in the heap.
     * @param node The handle of the node to update. Must be a valid handle obtained from insert().
     * @param from_value The original value of the node (used for offset calculation).
     * @param to_value The new, smaller value for the node. Must be less than or equal to from_value.
     */
    void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
        assert(node != nullptr && "Node handle cannot be null.");

        assert(to_value <= node->value && "New value must be smaller or equal to current node value.");

        ValueType additional_offset = from_value - node->value;
        node->child_offset += additional_offset;
        node->value = to_value;

        if (node == root_) {
            return;
        }

        if (node->left_up != nullptr) {

            Node* parent_or_left_sibling = node->left_up;
            if (parent_or_left_sibling->child == node) {

                parent_or_left_sibling->child = node->sibling;
            } else {

                parent_or_left_sibling->sibling = node->sibling;
            }

            if (node->sibling != nullptr) {
                node->sibling->left_up = parent_or_left_sibling;
            }

            node->left_up = nullptr;
            node->sibling = nullptr;

            root_ = link(root_, node);
            assert(root_ != nullptr);
        }

    }

    /**
     * @brief Deletes and returns the minimum element (root) from the heap.
     * @param[out] value Pointer to store the minimum value.
     * @param[out] payload Pointer to store the payload of the minimum element.
     * @return True if the heap was not empty and an element was deleted, false otherwise.
     */
    bool delete_min(ValueType* value, PayloadType* payload) {
        assert(value != nullptr && payload != nullptr);
        assert(buffer_ != nullptr && "Shared buffer required for delete_min.");
        if (root_ == nullptr) {
            return false;
        }

        Node* old_root = root_;
        *value = old_root->value;
        *payload = old_root->payload;

        buffer_->clear();
        Node* current_child = old_root->child;
        while (current_child != nullptr) {

            Node* next_sibling = current_child->sibling;

            current_child->value += old_root->child_offset;
            current_child->child_offset += old_root->child_offset;

            current_child->left_up = nullptr;
            current_child->sibling = nullptr;

            buffer_->push_back(current_child);
            current_child = next_sibling;
        }

        delete old_root;
        root_ = nullptr;

        if (buffer_->empty()) {
            return true;
        }

        size_t num_children = buffer_->size();
        size_t merged_children = 0;
        size_t write_idx = 0;
        while (merged_children + 1 < num_children) {
            (*buffer_)[write_idx] = link((*buffer_)[merged_children], (*buffer_)[merged_children + 1]);
            merged_children += 2;
            write_idx++;
        }

        if (merged_children < num_children) {
            (*buffer_)[write_idx] = (*buffer_)[merged_children];
            write_idx++;
        }

        buffer_->resize(write_idx);

        if (write_idx > 0) {
            root_ = (*buffer_)[write_idx - 1];
            for (int i = static_cast<int>(write_idx) - 2; i >= 0; --i) {
                root_ = link(root_, (*buffer_)[i]);
            }
        }

        assert((buffer_->empty() && root_ == nullptr) || (!buffer_->empty() && root_ != nullptr));

        buffer_->clear();
        return true;
    }

    /**
     * @brief Melds two pairing heaps into one.
     * The original heaps (heap1, heap2) become empty after the operation.
     * Note: Both heaps must have been constructed with the *same* shared buffer pointer.
     * @param heap1 The first heap to meld.
     * @param heap2 The second heap to meld.
     * @return A new PairingHeap containing elements from both heap1 and heap2.
     */
    [[nodiscard]] static PairingHeap meld(PairingHeap* heap1, PairingHeap* heap2) {
        assert(heap1 != nullptr && heap2 != nullptr);
        assert(heap1->buffer_ == heap2->buffer_ && "Heaps must share the same buffer for melding.");

        PairingHeap result(heap1->buffer_);
        result.root_ = link(heap1->root_, heap2->root_);

        heap1->root_ = nullptr;
        heap2->root_ = nullptr;

        return result;
    }

  private:
    Node* root_;
    std::vector<ItemHandle>* buffer_;

    /**
     * @brief Links two heap trees, maintaining the heap property.
     * The node with the smaller value becomes the root of the linked tree.
     * Adjusts values and offsets based on the linking process.
     * @param node1 Root of the first tree. Can be nullptr.
     * @param node2 Root of the second tree. Can be nullptr.
     * @return The root of the merged heap tree.
     */
    [[nodiscard]] static Node* link(Node* node1, Node* node2) {
        if (node1 == nullptr) return node2;
        if (node2 == nullptr) return node1;

        Node* smaller_node = node1;
        Node* larger_node = node2;

        if (node2->value < node1->value) {
            std::swap(smaller_node, larger_node);
        }
        assert(smaller_node->value <= larger_node->value);

        larger_node->sibling = smaller_node->child;
        if (smaller_node->child != nullptr) {
            assert(smaller_node->child->left_up == smaller_node);
            smaller_node->child->left_up = larger_node;
        }
        larger_node->left_up = smaller_node;
        smaller_node->child = larger_node;

        larger_node->value -= smaller_node->child_offset;

        larger_node->child_offset -= smaller_node->child_offset;

        return smaller_node;
    }
};

}