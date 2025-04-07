#pragma once

#include <string_view>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <limits>
#include <optional>
#include <functional>
#include <memory>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "logger.h"
#include "pairing_heap.h"
#include "priority_queue.h"

namespace cluster_approx {
    namespace internal {
        class IPruner;
        struct PruningContext;
    }

    class PCSTFast {
        public:
            using ValueType = double;
            using IndexType = int;
            using PayloadType = int;

            enum class PruningMethod {
                kNoPruning = 0,
                kSimplePruning,
                kGWPruning,
                kStrongPruning,
                kConnectFinalComponents,
                kUnknownPruning,
            };

            struct Statistics {
                long long total_num_edge_events = 0;
                long long num_deleted_edge_events = 0;
                long long num_merged_edge_events = 0;
                long long total_num_merge_events = 0;
                long long num_active_active_merge_events = 0;
                long long num_active_inactive_merge_events = 0;
                long long total_num_edge_growth_events = 0;
                long long num_active_active_edge_growth_events = 0;
                long long num_active_inactive_edge_growth_events = 0;
                long long num_cluster_events = 0;

                void reset() {
                    total_num_edge_events = 0;
                    num_deleted_edge_events = 0;
                    num_merged_edge_events = 0;
                    total_num_merge_events = 0;
                    num_active_active_merge_events = 0;
                    num_active_inactive_merge_events = 0;
                    total_num_edge_growth_events = 0;
                    num_active_active_edge_growth_events = 0;
                    num_active_inactive_edge_growth_events = 0;
                    num_cluster_events = 0;
                }
            };

            struct EdgeInfo {
                IndexType inactive_merge_event = kInvalidIndex;
            };

            struct InactiveMergeEvent {
                IndexType active_cluster_index = kInvalidIndex;
                IndexType inactive_cluster_index = kInvalidIndex;
                IndexType active_cluster_node = kInvalidIndex;
                IndexType inactive_cluster_node = kInvalidIndex;
            };

            struct Cluster {
                PairingHeap<ValueType, PayloadType> edge_parts;
                bool active = false;
                ValueType active_start_time = 0.0;
                ValueType active_end_time = -1.0;
                IndexType merged_into = kInvalidIndex;
                ValueType prize_sum = 0.0;
                ValueType subcluster_moat_sum = 0.0;
                ValueType moat = 0.0;
                bool contains_root = false;
                IndexType skip_up = kInvalidIndex;
                ValueType skip_up_sum = 0.0;
                IndexType merged_along = kInvalidIndex;
                IndexType child_cluster_1 = kInvalidIndex;
                IndexType child_cluster_2 = kInvalidIndex;
                bool necessary = false;

                Cluster() = default;
            };

            static constexpr IndexType kNoRoot = -1;
            static constexpr IndexType kInvalidIndex = -1;
            static constexpr ValueType kEpsilon = 1e-9;

        private:
            using PairingHeapType = PairingHeap<ValueType, PayloadType>;
            using PriorityQueueType = PriorityQueue<ValueType, IndexType>;

            struct EdgePart {
                ValueType next_event_val = 0.0;
                bool deleted = false;
                PairingHeapType::ItemHandle heap_node = PairingHeapType::kInvalidHandle;
            };

            const std::vector<std::pair<IndexType, IndexType>>& edges_;
            const std::vector<ValueType>& prizes_;
            const std::vector<ValueType>& costs_;
            const IndexType root_;
            const int target_num_active_clusters_;
            Statistics stats_;
            std::vector<EdgePart> edge_parts_;
            std::vector<EdgeInfo> edge_info_;
            std::vector<Cluster> clusters_;
            std::vector<InactiveMergeEvent> inactive_merge_events_;
            PriorityQueueType clusters_deactivation_;
            PriorityQueueType clusters_next_edge_event_;
            ValueType current_time_ = 0.0;
            std::unique_ptr<internal::IPruner> pruner_;
            internal::Logger& logger_;
            std::vector<uint8_t> node_good_;
            std::vector<int> phase1_result_;

            void initialize_clusters();
            void initialize_edges();

            struct NextEvent {
                ValueType time = std::numeric_limits<ValueType>::infinity();
                enum class Type { EDGE, CLUSTER_DEACTIVATION, NONE } type = Type::NONE;
                IndexType cluster_index = kInvalidIndex;
                PayloadType edge_part_index = kInvalidIndex;
            };

            [[nodiscard]] NextEvent get_next_event() const;
            void process_event(const NextEvent& event);
            void handle_edge_event(PayloadType edge_part_index);
            void handle_cluster_deactivation_event(IndexType cluster_index);

            struct EdgeProcessingInfo {
                ValueType sum_current = 0.0;
                ValueType sum_other = 0.0;
                ValueType finished_moat_current = 0.0;
                ValueType finished_moat_other = 0.0;
                IndexType current_cluster_idx = kInvalidIndex;
                IndexType other_cluster_idx = kInvalidIndex;
                IndexType edge_idx = kInvalidIndex;
                PayloadType current_part_idx = kInvalidIndex;
                PayloadType other_part_idx = kInvalidIndex;
                ValueType cost = 0.0;
                ValueType remainder = 0.0;
            };

            [[nodiscard]] std::optional<EdgeProcessingInfo> get_edge_processing_info(PayloadType edge_part_index);
            void merge_clusters(const EdgeProcessingInfo& info);
            void handle_active_active_growth(const EdgeProcessingInfo& info);
            void handle_active_inactive_growth(const EdgeProcessingInfo& info);
            void update_cluster_queues_post_merge(IndexType new_cluster_index, IndexType child1_index, IndexType child2_index);
            void update_cluster_queues_post_growth(IndexType cluster_idx);
            void log_inactive_merge_event(const EdgeProcessingInfo& info);
            void get_sum_on_edge_part(PayloadType edge_part_index, ValueType* total_sum,
                                    ValueType* finished_moat_sum, IndexType* current_cluster_index);
            [[nodiscard]] IndexType find_representative_and_compress(IndexType start_node, ValueType& path_sum_out);
            void select_initial_active_clusters(std::vector<uint8_t>& node_good_flag);

            [[nodiscard]] static constexpr PayloadType get_other_edge_part_index(PayloadType edge_part_index) {
                return edge_part_index ^ 1;
            }
        public:
            [[nodiscard]] static PruningMethod parse_pruning_method(std::string_view input);

            PCSTFast(const std::vector<std::pair<IndexType, IndexType>>& edges,
                    const std::vector<ValueType>& prizes,
                    const std::vector<ValueType>& costs,
                    IndexType root,
                    int target_num_active_clusters,
                    PruningMethod pruning,
                    internal::Logger& logger
                    );
            ~PCSTFast();
            PCSTFast(const PCSTFast&) = delete;
            PCSTFast& operator=(const PCSTFast&) = delete;
            PCSTFast(PCSTFast&&) = delete;
            PCSTFast& operator=(PCSTFast&&) = delete;
            [[nodiscard]] bool run(std::vector<IndexType>* result_nodes, std::vector<IndexType>* result_edges);
            void get_statistics(Statistics* s) const;
    };

    namespace internal {
         std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method);
    }

}
