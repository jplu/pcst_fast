#pragma once

#include <vector>
#include <deque>
#include <utility>
#include <limits>
#include <cassert>
#include <memory>
#include <span>

namespace cluster_approx {

/**
 * @brief Allocator for Pairing Heap nodes.
 * Uses std::deque to ensure pointer stability upon insertion (allocation),
 * preventing invalidation of existing node pointers when the container grows.
 */
template <typename ValueType, typename PayloadType>
class PairingHeapAllocator {
public:
    struct Node {
        Node* sibling = nullptr;
        Node* child = nullptr;
        Node* left_up = nullptr;
        ValueType value = ValueType{};
        ValueType child_offset = ValueType{};
        PayloadType payload = PayloadType{};
    };

    // Capacity is accepted for API compatibility but std::deque manages dynamic growth automatically.
    explicit PairingHeapAllocator(size_t /*capacity*/) {}

    Node* allocate(ValueType value, PayloadType payload) {
        nodes_.emplace_back();
        Node* node = &nodes_.back();
        node->value = value;
        node->payload = payload;
        return node;
    }

    // No deallocate needed for individual nodes during the run; 
    // memory is reclaimed when the allocator is destroyed.
    
    PairingHeapAllocator(PairingHeapAllocator&&) = default;
    PairingHeapAllocator& operator=(PairingHeapAllocator&&) = default;
    PairingHeapAllocator(const PairingHeapAllocator&) = delete;
    PairingHeapAllocator& operator=(const PairingHeapAllocator&) = delete;

private:
    std::deque<Node> nodes_;
};


/**
 * @brief Implements a Pairing Heap data structure.
 *
 * Uses an external allocator for node management.
 *
 * @tparam ValueType The type of the values (keys).
 * @tparam PayloadType The type of the payload.
 */
template <typename ValueType, typename PayloadType>
class PairingHeap {
  public:
    using AllocatorType = PairingHeapAllocator<ValueType, PayloadType>;
    using Node = typename AllocatorType::Node;
    using ItemHandle = Node*;

    /**
     * @brief Constructs a PairingHeap.
     * @param allocator Pointer to the shared allocator. Must outlive the heap.
     * @param shared_buffer Pointer to shared workspace buffer.
     */
    PairingHeap(AllocatorType* allocator, std::vector<ItemHandle>* shared_buffer) 
        : root_(nullptr), allocator_(allocator), buffer_(shared_buffer) {
        assert(allocator_ != nullptr && "Allocator cannot be null.");
        assert(buffer_ != nullptr && "Shared buffer cannot be null.");
    }

    PairingHeap(PairingHeap&& other) noexcept
        : root_(other.root_), allocator_(other.allocator_), buffer_(other.buffer_) {
        other.root_ = nullptr;
    }

    PairingHeap& operator=(PairingHeap&& other) noexcept {
        if (this != &other) {
            root_ = other.root_;
            allocator_ = other.allocator_;
            buffer_ = other.buffer_;
            other.root_ = nullptr;
        }
        return *this;
    }

    // Copying deleted
    PairingHeap(const PairingHeap&) = delete;
    PairingHeap& operator=(const PairingHeap&) = delete;

    [[nodiscard]] bool is_empty() const noexcept {
        return root_ == nullptr;
    }

    [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
        if (root_ != nullptr) {
            *value = root_->value;
            *payload = root_->payload;
            return true;
        }
        return false;
    }

    [[nodiscard]] ItemHandle insert(ValueType value, PayloadType payload) {
        Node* new_node = allocator_->allocate(value, payload);
        root_ = link(root_, new_node);
        return new_node;
    }

    void add_to_heap(ValueType value) {
        if (root_ != nullptr) {
            root_->value += value;
            root_->child_offset += value;
        }
    }

    void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
        assert(node != nullptr);
        assert(to_value <= node->value);

        ValueType additional_offset = from_value - node->value;
        node->child_offset += additional_offset;
        node->value = to_value;

        if (node == root_) return;

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
        }
    }

    bool delete_min(ValueType* value, PayloadType* payload) {
        if (root_ == nullptr) return false;

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

        // We do not delete old_root here; the allocator owns it.
        root_ = nullptr;

        if (buffer_->empty()) return true;

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
        buffer_->clear();
        return true;
    }

    [[nodiscard]] static PairingHeap meld(PairingHeap* heap1, PairingHeap* heap2) {
        assert(heap1->allocator_ == heap2->allocator_);
        assert(heap1->buffer_ == heap2->buffer_);

        PairingHeap result(heap1->allocator_, heap1->buffer_);
        result.root_ = link(heap1->root_, heap2->root_);

        heap1->root_ = nullptr;
        heap2->root_ = nullptr;

        return result;
    }

  private:
    Node* root_;
    AllocatorType* allocator_;
    std::vector<ItemHandle>* buffer_;

    static Node* link(Node* node1, Node* node2) {
        if (node1 == nullptr) return node2;
        if (node2 == nullptr) return node1;

        Node* smaller_node = node1;
        Node* larger_node = node2;

        if (node2->value < node1->value) {
            std::swap(smaller_node, larger_node);
        }

        larger_node->sibling = smaller_node->child;
        if (smaller_node->child != nullptr) {
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