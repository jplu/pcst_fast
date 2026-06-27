#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <cassert>

namespace cluster_approx {

/**
 * @brief Allocator for Pairing Heap nodes.
 * Uses a contiguous std::vector to guarantee cache-friendly contiguous layout.
 * Since indices (int32_t) are used instead of raw pointers, reallocations during
 * vector growth do not invalidate structural heap handles.
 */
template <typename ValueType, typename PayloadType>
class PairingHeapAllocator {
public:
    struct Node {
        int32_t sibling = -1;
        int32_t child = -1;
        int32_t left_up = -1;
        ValueType value = ValueType{};
        ValueType child_offset = ValueType{};
        PayloadType payload = PayloadType{};
    };

    explicit PairingHeapAllocator(size_t capacity) {
        nodes_.reserve(capacity);
    }

    int32_t allocate(ValueType value, PayloadType payload) {
        int32_t idx = static_cast<int32_t>(nodes_.size());
        nodes_.push_back(Node{
            .sibling = -1,
            .child = -1,
            .left_up = -1,
            .value = value,
            .child_offset = ValueType{},
            .payload = payload
        });
        return idx;
    }

    Node& operator[](int32_t idx) noexcept {
        return nodes_[idx];
    }

    const Node& operator[](int32_t idx) const noexcept {
        return nodes_[idx];
    }

    size_t size() const noexcept { return nodes_.size(); }

    void clear() noexcept { nodes_.clear(); }

private:
    std::vector<Node> nodes_;
};


/**
 * @brief Implements an index-based Pairing Heap data structure.
 *
 * Uses an external contiguous allocator for node management.
 *
 * @tparam ValueType The type of the values (keys).
 * @tparam PayloadType The type of the payload.
 */
template <typename ValueType, typename PayloadType>
class PairingHeap {
public:
    using AllocatorType = PairingHeapAllocator<ValueType, PayloadType>;
    using Node = typename AllocatorType::Node;
    using ItemHandle = int32_t;

    /**
     * @brief Constructs a PairingHeap.
     * @param allocator Pointer to the shared allocator. Must outlive the heap.
     * @param shared_buffer Pointer to shared workspace buffer.
     */
    PairingHeap(AllocatorType* allocator, std::vector<ItemHandle>* shared_buffer) 
        : root_(-1), allocator_(allocator), buffer_(shared_buffer) {
        assert(allocator_ != nullptr && "Allocator cannot be null.");
        assert(buffer_ != nullptr && "Shared buffer cannot be null.");
    }

    PairingHeap(PairingHeap&& other) noexcept
        : root_(other.root_), allocator_(other.allocator_), buffer_(other.buffer_) {
        other.root_ = -1;
    }

    PairingHeap& operator=(PairingHeap&& other) noexcept {
        if (this != &other) {
            root_ = other.root_;
            allocator_ = other.allocator_;
            buffer_ = other.buffer_;
            other.root_ = -1;
        }
        return *this;
    }

    PairingHeap(const PairingHeap&) = delete;
    PairingHeap& operator=(const PairingHeap&) = delete;

    [[nodiscard]] bool is_empty() const noexcept {
        return root_ == -1;
    }

    [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
        if (root_ != -1) {
            const auto& root_node = (*allocator_)[root_];
            *value = root_node.value;
            *payload = root_node.payload;
            return true;
        }
        return false;
    }

    [[nodiscard]] ItemHandle insert(ValueType value, PayloadType payload) {
        int32_t new_node = allocator_->allocate(value, payload);
        root_ = link(root_, new_node);
        return new_node;
    }

    void add_to_heap(ValueType value) {
        if (root_ != -1) {
            auto& root_node = (*allocator_)[root_];
            root_node.value += value;
            root_node.child_offset += value;
        }
    }

    void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
        assert(node != -1);
        auto& n = (*allocator_)[node];
        assert(to_value <= n.value);

        ValueType additional_offset = from_value - n.value;
        n.child_offset += additional_offset;
        n.value = to_value;

        if (node == root_) return;

        if (n.left_up != -1) {
            int32_t parent_or_left_sibling = n.left_up;
            auto& p_or_l = (*allocator_)[parent_or_left_sibling];
            if (p_or_l.child == node) {
                p_or_l.child = n.sibling;
            } else {
                p_or_l.sibling = n.sibling;
            }

            if (n.sibling != -1) {
                (*allocator_)[n.sibling].left_up = parent_or_left_sibling;
            }

            n.left_up = -1;
            n.sibling = -1;

            root_ = link(root_, node);
        }
    }

    bool delete_min(ValueType* value, PayloadType* payload) {
        if (root_ == -1) return false;

        auto& old_root = (*allocator_)[root_];
        *value = old_root.value;
        *payload = old_root.payload;

        buffer_->clear();
        int32_t current_child = old_root.child;
        ValueType root_offset = old_root.child_offset;
        
        while (current_child != -1) {
            auto& child_node = (*allocator_)[current_child];
            int32_t next_sibling = child_node.sibling;
            child_node.value += root_offset;
            child_node.child_offset += root_offset;
            child_node.left_up = -1;
            child_node.sibling = -1;
            buffer_->push_back(current_child);
            current_child = next_sibling;
        }

        root_ = -1;

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
        result.root_ = link(heap1->root_, heap2->root_, heap1->allocator_);

        heap1->root_ = -1;
        heap2->root_ = -1;

        return result;
    }

private:
    int32_t root_;
    AllocatorType* allocator_;
    std::vector<ItemHandle>* buffer_;

    static int32_t link(int32_t node1, int32_t node2, AllocatorType* allocator) noexcept {
        if (node1 == -1) return node2;
        if (node2 == -1) return node1;

        int32_t smaller_node = node1;
        int32_t larger_node = node2;

        auto& n1 = (*allocator)[node1];
        auto& n2 = (*allocator)[node2];

        if (n2.value < n1.value) {
            std::swap(smaller_node, larger_node);
        }

        auto& s_node = (*allocator)[smaller_node];
        auto& l_node = (*allocator)[larger_node];

        l_node.sibling = s_node.child;
        if (s_node.child != -1) {
            (*allocator)[s_node.child].left_up = larger_node;
        }
        l_node.left_up = smaller_node;
        s_node.child = larger_node;

        l_node.value -= s_node.child_offset;
        l_node.child_offset -= s_node.child_offset;

        return smaller_node;
    }

    int32_t link(int32_t node1, int32_t node2) noexcept {
        return link(node1, node2, allocator_);
    }
};

}