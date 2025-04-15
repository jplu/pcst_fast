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
        logger_->log(LogLevel::WARNING, "Attempted to mark necessary from invalid node/cluster index: {}", start_node_index);
        return;
    }

    logger_->log(LogLevel::TRACE, "  Attempting to mark necessary chain starting from original node index: {}", start_node_index);

    while (static_cast<size_t>(current_cluster_idx) < clusters.size() &&
            !clusters[current_cluster_idx].necessary) {
        logger_->log(LogLevel::TRACE, "    Marking cluster {} as necessary.", current_cluster_idx);
        clusters[current_cluster_idx].necessary = true;
        if (clusters[current_cluster_idx].merged_into != kInvalidClusterId) {
            current_cluster_idx = clusters[current_cluster_idx].merged_into;
        } else {

            logger_->log(LogLevel::TRACE, "    Reached merge tree root (cluster {}), stopping necessary propagation.", current_cluster_idx);
            return;
        }
    }

    if(static_cast<size_t>(current_cluster_idx) < clusters.size() && clusters[current_cluster_idx].necessary) {
        logger_->log(LogLevel::TRACE, "    Stopped necessary propagation at cluster {} (already marked).", current_cluster_idx);
    }
}

PruningResult GWPruner::prune(const PruningInput& input) {
    input_ = &input;
    logger_ = input.logger;
    assert(logger_ != nullptr);
    num_nodes_ = input.graph.prizes.size();
    node_deleted_.assign(num_nodes_, false);
    node_queue_.clear();

    clusters_ptr_ = const_cast<std::vector<Cluster>*>(&input.core_result.final_cluster_state);
    assert(clusters_ptr_ != nullptr);

    logger_->log(LogLevel::INFO, "Applying GWPruning strategy.");

    std::vector<EdgeId> intermediate_edges;
    intermediate_edges.reserve(input.core_result.phase1_edges.size());
    for (EdgeId edge_idx : input.core_result.phase1_edges) {
        assert(static_cast<size_t>(edge_idx) < input.graph.edges.size());
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;
        assert(static_cast<size_t>(u) < num_nodes_ && static_cast<size_t>(v) < num_nodes_);
        assert(static_cast<size_t>(u) < input.core_result.initial_node_filter.size());
        assert(static_cast<size_t>(v) < input.core_result.initial_node_filter.size());

        if (input.core_result.initial_node_filter[u] && input.core_result.initial_node_filter[v]) {
            intermediate_edges.push_back(edge_idx);
        }
    }
    logger_->log(LogLevel::DEBUG, "GWPruning: Starting with {} intermediate edges (filtered from phase1).", intermediate_edges.size());

    if (intermediate_edges.empty()) {
        logger_->log(LogLevel::INFO, "No intermediate edges after filtering, GW pruning results in empty graph.");

        return { build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter), {} };
    }

    neighbors_ = build_adjacency_list(intermediate_edges, input.graph);
    logger_->log(LogLevel::DEBUG, "Built adjacency list for GW pruning graph.");

    std::vector<EdgeId> final_edges;
    final_edges.reserve(intermediate_edges.size());

    for (int i = static_cast<int>(intermediate_edges.size()) - 1; i >= 0; --i) {
        EdgeId current_edge_index = intermediate_edges[i];
        NodeId u = input.graph.edges[current_edge_index].first;
        NodeId v = input.graph.edges[current_edge_index].second;

        logger_->log(LogLevel::TRACE, "Processing edge {} ({}, {}) in reverse order.", current_edge_index, u, v);

        if (node_deleted_[u] && node_deleted_[v]) {
            logger_->log(LogLevel::TRACE, "  Skipping edge {} ({},{}): Both endpoints already deleted.", current_edge_index, u, v);
            continue;
        }

        EventId merge_event_id = input.core_result.edge_inactive_merge_event_ids[current_edge_index];

        if (merge_event_id == kInvalidEventId) {

            logger_->log(LogLevel::TRACE, "  Edge {} ({},{}) from Active-Active merge. Keeping.", current_edge_index, u, v);
            final_edges.push_back(current_edge_index);

            mark_clusters_as_necessary_from_node(u);
            mark_clusters_as_necessary_from_node(v);
        } else {

            logger_->log(LogLevel::TRACE, "  Edge {} ({},{}) from Active-Inactive merge (EventID: {}).", current_edge_index, u, v, merge_event_id);
            assert(static_cast<size_t>(merge_event_id) < input.core_result.inactive_merge_events.size());
            const InactiveMergeEvent& merge_event = input.core_result.inactive_merge_events[merge_event_id];

            NodeId active_side_node = merge_event.active_cluster_node;
            NodeId inactive_side_node = merge_event.inactive_cluster_node;
            ClusterId inactive_cluster_index = merge_event.inactive_cluster_index;

            assert(static_cast<size_t>(inactive_cluster_index) < clusters_ptr_->size() && "Inactive cluster index out of bounds!");

            if ((*clusters_ptr_)[inactive_cluster_index].necessary) {
                logger_->log(LogLevel::TRACE, "  Inactive cluster index {} is necessary. Keeping edge {}.", inactive_cluster_index, current_edge_index);
                final_edges.push_back(current_edge_index);

                mark_clusters_as_necessary_from_node(active_side_node);
                mark_clusters_as_necessary_from_node(inactive_side_node);
            } else {

                bool inactive_is_root = (inactive_side_node == input_->graph.root && input_->graph.root != kInvalidNodeId);
                if (inactive_is_root) {
                    logger_->log(LogLevel::TRACE, "  Inactive side node {} is root. Keeping edge {}.", inactive_side_node, current_edge_index);
                    final_edges.push_back(current_edge_index);
                    mark_clusters_as_necessary_from_node(active_side_node);
                    mark_clusters_as_necessary_from_node(inactive_side_node);
                } else {
                    logger_->log(LogLevel::TRACE, "  Inactive cluster index {} (node {}) is not necessary/root. Discarding edge {} and marking inactive side node {} deleted.",
                                 inactive_cluster_index, inactive_side_node, current_edge_index, inactive_side_node);
                    mark_nodes_as_deleted(inactive_side_node, active_side_node);
                }
            }
        }
    }

    std::reverse(final_edges.begin(), final_edges.end());
    logger_->log(LogLevel::DEBUG, "GWPruning: Selected {} final edges.", final_edges.size());

    PruningResult result;
    result.edges = std::move(final_edges);

    result.nodes = build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter);

    logger_->log(LogLevel::DEBUG, "GWPruning: Derived {} final nodes.", result.nodes.size());
    logger_->log(LogLevel::INFO, "GWPruning completed.");

    return result;
}

void GWPruner::mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index) {
    logger_->log(LogLevel::TRACE, "  Marking deleted starting from node {}, parent {}", start_node_index, parent_node_index);
    assert(static_cast<size_t>(start_node_index) < num_nodes_);
    assert(parent_node_index == kInvalidNodeId || static_cast<size_t>(parent_node_index) < num_nodes_);

    node_queue_.clear();

    if (node_deleted_[start_node_index]) {
        logger_->log(LogLevel::TRACE, "    Node {} was already marked deleted.", start_node_index);
        return;
    }

    node_deleted_[start_node_index] = true;
    node_queue_.push_back(start_node_index);
    logger_->log(LogLevel::TRACE, "    Marked node {} as deleted.", start_node_index);

    size_t current_idx = 0;
    while(current_idx < node_queue_.size()) {
        NodeId current_node = node_queue_[current_idx++];

        if (static_cast<size_t>(current_node) >= neighbors_.size()) {
            logger_->log(LogLevel::WARNING, "    Node {} out of bounds for neighbors_ list in mark_nodes_as_deleted.", current_node);
            continue;
        }

        for(const auto& edge_pair : neighbors_[current_node]) {
            NodeId neighbor_node = edge_pair.first;

            if (static_cast<size_t>(neighbor_node) >= node_deleted_.size()) {
                logger_->log(LogLevel::WARNING, "    Neighbor node {} out of bounds for node_deleted_ check in mark_nodes_as_deleted.", neighbor_node);
                continue;
            }

            if (neighbor_node == parent_node_index) {
                continue;
            }

            if (!node_deleted_[neighbor_node]) {
                node_deleted_[neighbor_node] = true;
                node_queue_.push_back(neighbor_node);
                logger_->log(LogLevel::TRACE, "    Marked node {} as deleted (neighbor of {}).", neighbor_node, current_node);
            }
        }

    }
}

}
}