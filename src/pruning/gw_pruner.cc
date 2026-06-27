#include "pcst_fast/pruning/gw_pruner.h"
#include "pcst_fast/pruning/pruning_utils.h"
#include "pcst_fast/pcst_core_internals.h"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>

namespace cluster_approx {

using cluster_approx::Cluster;

namespace pruning {

void GWPruner::mark_clusters_as_necessary_from_node(NodeId start_node_index) {
    assert(clusters_ptr_ != nullptr && "Cluster pointer must be set before marking necessary.");
    auto& clusters = *clusters_ptr_;

    ClusterId current_cluster_idx = start_node_index;

    if (current_cluster_idx < 0 || static_cast<size_t>(current_cluster_idx) >= clusters.size()) {
        return;
    }

    while (static_cast<size_t>(current_cluster_idx) < clusters.size() &&
            !clusters[current_cluster_idx].necessary) {
        clusters[current_cluster_idx].necessary = true;
        if (clusters[current_cluster_idx].merged_into != kInvalidClusterId) {
            current_cluster_idx = clusters[current_cluster_idx].merged_into;
        } else {
            return;
        }
    }
}

PruningResult GWPruner::prune(const PruningInput& input) {
    input_ = &input;
    logger_ = input.logger;
    assert(logger_ != nullptr);
    num_nodes_ = input.graph.prizes.size();
    node_deleted_.assign(num_nodes_, 0);
    node_queue_.clear();

    clusters_ptr_ = &input.core_result.final_cluster_state;
    assert(clusters_ptr_ != nullptr);

    logger_->log(LogLevel::INFO, "Applying GWPruning strategy.");

    std::vector<EdgeId> intermediate_edges;
    intermediate_edges.reserve(input.core_result.phase1_edges.size());
    for (EdgeId edge_idx : input.core_result.phase1_edges) {
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;

        if (input.core_result.initial_node_filter[u] && input.core_result.initial_node_filter[v]) {
            intermediate_edges.push_back(edge_idx);
        }
    }

    if (intermediate_edges.empty()) {
        return { build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter), {} };
    }

    neighbors_ = build_adjacency_list_csr(intermediate_edges, input.graph);

    std::vector<EdgeId> final_edges;
    final_edges.reserve(intermediate_edges.size());

    for (int i = static_cast<int>(intermediate_edges.size()) - 1; i >= 0; --i) {
        EdgeId current_edge_index = intermediate_edges[i];
        NodeId u = input.graph.edges[current_edge_index].first;
        NodeId v = input.graph.edges[current_edge_index].second;

        if (node_deleted_[u] && node_deleted_[v]) {
            continue;
        }

        EventId merge_event_id = input.core_result.edge_inactive_merge_event_ids[current_edge_index];

        if (merge_event_id == kInvalidEventId) {
            final_edges.push_back(current_edge_index);
            mark_clusters_as_necessary_from_node(u);
            mark_clusters_as_necessary_from_node(v);
        } else {
            const InactiveMergeEvent& merge_event = input.core_result.inactive_merge_events[merge_event_id];
            NodeId active_side_node = merge_event.active_cluster_node;
            NodeId inactive_side_node = merge_event.inactive_cluster_node;
            ClusterId inactive_cluster_index = merge_event.inactive_cluster_index;

            if ((*clusters_ptr_)[inactive_cluster_index].necessary) {
                final_edges.push_back(current_edge_index);
                mark_clusters_as_necessary_from_node(active_side_node);
                mark_clusters_as_necessary_from_node(inactive_side_node);
            } else {
                bool inactive_is_root = (inactive_side_node == input_->graph.root && input_->graph.root != kInvalidNodeId);
                if (inactive_is_root) {
                    final_edges.push_back(current_edge_index);
                    mark_clusters_as_necessary_from_node(active_side_node);
                    mark_clusters_as_necessary_from_node(inactive_side_node);
                } else {
                    mark_nodes_as_deleted(inactive_side_node, active_side_node);
                }
            }
        }
    }

    std::reverse(final_edges.begin(), final_edges.end());

    PruningResult result;
    result.edges = std::move(final_edges);
    result.nodes = build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter);

    return result;
}

void GWPruner::mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index) {
    node_queue_.clear();

    if (node_deleted_[start_node_index]) {
        return;
    }

    node_deleted_[start_node_index] = 1;
    node_queue_.push_back(start_node_index);

    size_t current_idx = 0;
    while(current_idx < node_queue_.size()) {
        NodeId current_node = node_queue_[current_idx++];

        for(const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
            NodeId neighbor_node = edge_pair.first;

            if (neighbor_node == parent_node_index) {
                continue;
            }

            if (!node_deleted_[neighbor_node]) {
                node_deleted_[neighbor_node] = 1;
                node_queue_.push_back(neighbor_node);
            }
        }
    }
}

}
}