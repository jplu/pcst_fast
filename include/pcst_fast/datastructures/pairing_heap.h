/*#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <cassert>
#include <memory>
#include "pcst_fast/pcst_types.h"

namespace cluster_approx {

template <typename ValueType, typename PayloadType>
class PairingHeapAllocator {
public:
    struct alignas(32) Node {
        ValueType value;
        ValueType child_offset;
        int32_t sibling;
        int32_t child;
        int32_t left_up;
        PayloadType payload;
    };

private:
    static constexpr size_t ChunkShift = 16;
    static constexpr size_t ChunkSize = 1ULL << ChunkShift;
    static constexpr size_t ChunkMask = ChunkSize - 1;

    std::vector<std::unique_ptr<Node[]>> chunks_;
    size_t size_ = 0;

public:
    explicit PairingHeapAllocator(size_t capacity = 0) {
        if (capacity > 0) {
            size_t num_chunks = (capacity + ChunkSize - 1) >> ChunkShift;
            chunks_.reserve(num_chunks);
            for (size_t i = 0; i < num_chunks; ++i) {
                chunks_.push_back(std::make_unique<Node[]>(ChunkSize));
            }
        }
    }

    FORCE_INLINE int32_t allocate(ValueType value, PayloadType payload) {
        size_t idx = size_++;
        size_t chunk_idx = idx >> ChunkShift;
        if (chunk_idx >= chunks_.size()) {
            chunks_.push_back(std::make_unique<Node[]>(ChunkSize));
        }
        Node& node = chunks_[chunk_idx][idx & ChunkMask];
        node.value = value;
        node.child_offset = ValueType{};
        node.sibling = -1;
        node.child = -1;
        node.left_up = -1;
        node.payload = payload;
        return static_cast<int32_t>(idx);
    }

    FORCE_INLINE Node& operator[](int32_t idx) noexcept { 
        return chunks_[idx >> ChunkShift][idx & ChunkMask]; 
    }
    FORCE_INLINE const Node& operator[](int32_t idx) const noexcept { 
        return chunks_[idx >> ChunkShift][idx & ChunkMask]; 
    }

    size_t size() const noexcept { return size_; }
    void clear() noexcept { size_ = 0; } // Drops bounds lazily without free payload overhead
};

template <typename ValueType, typename PayloadType>
class PairingHeap {
public:
    using AllocatorType = PairingHeapAllocator<ValueType, PayloadType>;
    using Node = typename AllocatorType::Node;
    using ItemHandle = int32_t;

    PairingHeap(AllocatorType* allocator, std::vector<ItemHandle>* shared_buffer) 
        : root_(-1), allocator_(allocator), buffer_(shared_buffer) {}

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

    FORCE_INLINE [[nodiscard]] bool is_empty() const noexcept { return root_ == -1; }

    FORCE_INLINE [[nodiscard]] bool get_min(ValueType* value, PayloadType* payload) const {
        if (root_ != -1) {
            const auto& root_node = (*allocator_)[root_];
            *value = root_node.value;
            *payload = root_node.payload;
            return true;
        }
        return false;
    }

    FORCE_INLINE ItemHandle insert(ValueType value, PayloadType payload) {
        int32_t new_node = allocator_->allocate(value, payload);
        root_ = link(root_, new_node);
        return new_node;
    }

    FORCE_INLINE void add_to_heap(ValueType value) {
        if (root_ != -1 && value > 0.0) {
            auto& root_node = (*allocator_)[root_];
            root_node.value += value;
            root_node.child_offset += value;
        }
    }

    void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
        assert(node != -1);
        auto& n = (*allocator_)[node];

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

    FORCE_INLINE [[nodiscard]] static PairingHeap meld(PairingHeap* heap1, PairingHeap* heap2) {
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

    FORCE_INLINE static int32_t link(int32_t node1, int32_t node2, AllocatorType* allocator) noexcept {
        if (node1 == -1) return node2;
        if (node2 == -1) return node1;

        Node* s_node = &(*allocator)[node1];
        Node* l_node = &(*allocator)[node2];
        int32_t smaller_node = node1;
        int32_t larger_node = node2;

        if (l_node->value < s_node->value) {
            std::swap(smaller_node, larger_node);
            std::swap(s_node, l_node);
        }

        l_node->sibling = s_node->child;
        if (s_node->child != -1) {
            (*allocator)[s_node->child].left_up = larger_node;
        }
        l_node->left_up = smaller_node;
        s_node->child = larger_node;

        l_node->value -= s_node->child_offset;
        l_node->child_offset -= s_node->child_offset;

        return smaller_node;
    }

    FORCE_INLINE int32_t link(int32_t node1, int32_t node2) noexcept {
        return link(node1, node2, allocator_);
    }
};

}*/