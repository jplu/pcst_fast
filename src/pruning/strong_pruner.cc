#include "pcst_fast/pruning/strong_pruner.h"
#include "pcst_fast/pruning/pruning_utils.h"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>

namespace cluster_approx {
namespace pruning {

PruningResult StrongPruner::prune(const PruningInput& input) {

    input_ = &input;
    logger_ = input.logger;
    assert(logger_ != nullptr);
    num_nodes_ = input.graph.prizes.size();

    logger_->log(LogLevel::INFO, "Applying StrongPruning strategy.");

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
    logger_->log(LogLevel::DEBUG, "StrongPruning: Starting with {} intermediate edges (filtered from phase1).", intermediate_edges.size());

    if (intermediate_edges.empty()) {
        logger_->log(LogLevel::INFO, "No intermediate edges after filtering, Strong pruning results in empty graph.");

        std::vector<bool> temp_deleted(num_nodes_, false);
        return { build_final_node_set(num_nodes_, temp_deleted, input.core_result.initial_node_filter), {} };
    }

    node_deleted_.assign(num_nodes_, false);
    final_component_label_.assign(num_nodes_, kInvalidClusterId);
    final_components_.clear();
    root_component_index_ = kInvalidClusterId;
    strong_pruning_parent_.assign(num_nodes_, {kInvalidNodeId, 0.0});
    strong_pruning_payoff_.assign(num_nodes_, -1.0);
    dfs_stack_.clear();
    dfs_stack2_.clear();
    node_queue_.clear();

    neighbors_ = build_adjacency_list(intermediate_edges, input.graph);
    logger_->log(LogLevel::DEBUG, "Built adjacency list for Strong pruning graph.");

    logger_->log(LogLevel::DEBUG, "Finding connected components...");
    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes_); ++i) {

        bool node_is_relevant = !neighbors_[i].empty() ||
                                (static_cast<size_t>(i) < input.core_result.initial_node_filter.size() && input.core_result.initial_node_filter[i]);

        if (node_is_relevant && final_component_label_[i] == kInvalidClusterId) {
            final_components_.emplace_back();
            ClusterId current_component_idx = final_components_.size() - 1;
            logger_->log(LogLevel::TRACE, "Starting DFS for component {} from node {}.", current_component_idx, i);
            label_final_component(i, current_component_idx);
            logger_->log(LogLevel::DEBUG, "Found component {} with {} nodes.", current_component_idx, final_components_.back().size());
        }
    }
    logger_->log(LogLevel::INFO, "Identified {} connected components.", final_components_.size());

    for (ClusterId comp_idx = 0; comp_idx < static_cast<ClusterId>(final_components_.size()); ++comp_idx) {
        logger_->log(LogLevel::DEBUG, "Processing component {} (size {}).", comp_idx, final_components_[comp_idx].size());

        if (final_components_[comp_idx].empty()) {
            logger_->log(LogLevel::WARNING, "Component {} is empty, skipping.", comp_idx);
            continue;
        }

        if (comp_idx == root_component_index_) {

            logger_->log(LogLevel::DEBUG, "Component {} contains root {}. Pruning starting from root.", comp_idx, input_->graph.root);
            assert(input_->graph.root != kInvalidNodeId);
            strong_pruning_dfs(input_->graph.root, true);
        } else {

            logger_->log(LogLevel::DEBUG, "Component {} does not contain root. Finding best root...", comp_idx);
            NodeId best_root = find_best_component_root(comp_idx);
            logger_->log(LogLevel::DEBUG, "Best root for component {} is {}. Pruning from best root.", comp_idx, best_root);

            strong_pruning_parent_.assign(num_nodes_, {kInvalidNodeId, 0.0});
            strong_pruning_payoff_.assign(num_nodes_, -1.0);
            strong_pruning_dfs(best_root, true);
        }
    }

    std::vector<EdgeId> final_edges;
    final_edges.reserve(intermediate_edges.size());
    for (EdgeId edge_idx : intermediate_edges) {
        assert(static_cast<size_t>(edge_idx) < input.graph.edges.size());
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;
        assert(static_cast<size_t>(u) < num_nodes_ && static_cast<size_t>(v) < num_nodes_);

        if (!node_deleted_[u] && !node_deleted_[v]) {
            final_edges.push_back(edge_idx);
        } else {
            logger_->log(LogLevel::TRACE,"Discarding edge {} ({},{}) due to deleted endpoint(s).", edge_idx, u, v);
        }
    }
    logger_->log(LogLevel::DEBUG, "StrongPruning: Selected {} final edges.", final_edges.size());

    PruningResult result;
    result.edges = std::move(final_edges);

    result.nodes = build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter);

    logger_->log(LogLevel::DEBUG, "StrongPruning: Derived {} final nodes.", result.nodes.size());
    logger_->log(LogLevel::INFO, "StrongPruning completed.");

    input_ = nullptr;
    logger_ = nullptr;

    return result;
}

void StrongPruner::label_final_component(NodeId start_node_index, ClusterId component_index) {
    assert(static_cast<size_t>(start_node_index) < num_nodes_);
    assert(static_cast<size_t>(component_index) < final_components_.size());
    assert(final_component_label_[start_node_index] == kInvalidClusterId);

    dfs_stack2_.clear();
    dfs_stack2_.push_back(start_node_index);
    final_component_label_[start_node_index] = component_index;

    while (!dfs_stack2_.empty()) {
        NodeId current_node = dfs_stack2_.back();
        dfs_stack2_.pop_back();

        final_components_[component_index].push_back(current_node);
        if (current_node == input_->graph.root) {
            root_component_index_ = component_index;
            logger_->log(LogLevel::DEBUG, "  Node {} is root, marking component {} as root component.", current_node, component_index);
        }

        assert(static_cast<size_t>(current_node) < neighbors_.size());
        for (const auto& edge_pair : neighbors_[current_node]) {
            NodeId neighbor_node = edge_pair.first;
            assert(static_cast<size_t>(neighbor_node) < num_nodes_);
            if (final_component_label_[neighbor_node] == kInvalidClusterId) {
                final_component_label_[neighbor_node] = component_index;
                dfs_stack2_.push_back(neighbor_node);
            } else {

                assert(final_component_label_[neighbor_node] == component_index && "Node labelled with different component index during DFS!");
            }
        }
    }
}

void StrongPruner::strong_pruning_dfs(NodeId start_node_index, bool mark_as_deleted) {
    logger_->log(LogLevel::TRACE, "Entering strong_pruning_dfs from {}, mark_deleted={}", start_node_index, mark_as_deleted);
    assert(static_cast<size_t>(start_node_index) < num_nodes_);

    dfs_stack_.clear();

    strong_pruning_parent_[start_node_index] = {kInvalidNodeId, 0.0};
    dfs_stack_.push_back({true, start_node_index});

    while (!dfs_stack_.empty()) {
        auto [is_entry_call, current_node] = dfs_stack_.back();
        dfs_stack_.pop_back();

        if (is_entry_call) {
            logger_->log(LogLevel::TRACE, "  DFS Entry: Node {}", current_node);

            dfs_stack_.push_back({false, current_node});

            assert(static_cast<size_t>(current_node) < neighbors_.size());
            assert(static_cast<size_t>(current_node) < strong_pruning_parent_.size());

            for (const auto& edge_pair : neighbors_[current_node]) {
                NodeId neighbor_node = edge_pair.first;
                double edge_cost = edge_pair.second;

                if (neighbor_node == strong_pruning_parent_[current_node].first) {
                    continue;
                }

                logger_->log(LogLevel::TRACE, "    Pushing neighbor {} (cost {}) for entry.", neighbor_node, edge_cost);
                assert(static_cast<size_t>(neighbor_node) < strong_pruning_parent_.size());
                strong_pruning_parent_[neighbor_node] = {current_node, edge_cost};
                dfs_stack_.push_back({true, neighbor_node});
            }
        } else {
            logger_->log(LogLevel::TRACE, "  DFS Exit: Node {}", current_node);

            assert(static_cast<size_t>(current_node) < input_->graph.prizes.size());
            assert(static_cast<size_t>(current_node) < strong_pruning_payoff_.size());

            strong_pruning_payoff_[current_node] = input_->graph.prizes[current_node];
            logger_->log(LogLevel::TRACE, "    Node {} base prize: {:.4f}", current_node, strong_pruning_payoff_[current_node]);

            assert(static_cast<size_t>(current_node) < neighbors_.size());
            for (const auto& edge_pair : neighbors_[current_node]) {
                NodeId neighbor_node = edge_pair.first;
                double edge_cost = edge_pair.second;

                assert(static_cast<size_t>(neighbor_node) < strong_pruning_parent_.size());

                if (strong_pruning_parent_[neighbor_node].first != current_node) {
                    continue;
                }

                assert(static_cast<size_t>(neighbor_node) < strong_pruning_payoff_.size());
                assert(strong_pruning_payoff_[neighbor_node] >= -1.0 + 1e-9);

                double child_net_payoff = strong_pruning_payoff_[neighbor_node] - edge_cost;
                logger_->log(LogLevel::TRACE, "    Considering child {}: Payoff={:.4f}, EdgeCost={:.4f}, NetPayoff={:.4f}",
                             neighbor_node, strong_pruning_payoff_[neighbor_node], edge_cost, child_net_payoff);

                if (child_net_payoff <= 1e-9) {
                    if (mark_as_deleted) {
                        logger_->log(LogLevel::DEBUG, "    Pruning subtree at node {} (from parent {}): Net payoff {:.4f} <= 0.",
                                     neighbor_node, current_node, child_net_payoff);

                        mark_nodes_as_deleted(neighbor_node, current_node);
                    } else {
                        logger_->log(LogLevel::TRACE, "    Subtree at node {} would be pruned (payoff {:.4f}), but not marking.",
                                     neighbor_node, child_net_payoff);
                    }
                } else {

                    strong_pruning_payoff_[current_node] += child_net_payoff;
                    logger_->log(LogLevel::TRACE, "    Adding positive payoff {:.4f} from child {} to parent {}. New parent payoff: {:.4f}",
                                 child_net_payoff, neighbor_node, current_node, strong_pruning_payoff_[current_node]);
                }
            }
            logger_->log(LogLevel::TRACE, "  DFS Exit: Node {} final payoff: {:.4f}", current_node, strong_pruning_payoff_[current_node]);
        }
    }
    logger_->log(LogLevel::TRACE,"Exiting strong_pruning_dfs from {}", start_node_index);
}

NodeId StrongPruner::find_best_component_root(ClusterId component_index) {
    assert(static_cast<size_t>(component_index) < final_components_.size());
    const auto& component_nodes = final_components_[component_index];
    assert(!component_nodes.empty());

    NodeId initial_root = component_nodes[0];
    logger_->log(LogLevel::TRACE,"Finding best root for component {}, starting initial DFS from {}", component_index, initial_root);

    strong_pruning_parent_.assign(num_nodes_, {kInvalidNodeId, 0.0});
    strong_pruning_payoff_.assign(num_nodes_, -1.0);
    strong_pruning_dfs(initial_root, false);

    NodeId current_best_root = initial_root;
    double current_best_value = strong_pruning_payoff_[initial_root];
    logger_->log(LogLevel::TRACE,"  Initial root {} has payoff {:.4f}", initial_root, current_best_value);

    dfs_stack2_.clear();

    assert(static_cast<size_t>(initial_root) < neighbors_.size());
    for (const auto& edge_pair : neighbors_[initial_root]) {
        NodeId neighbor_node = edge_pair.first;

        assert(static_cast<size_t>(neighbor_node) < final_component_label_.size());
        if (final_component_label_[neighbor_node] == component_index) {
            dfs_stack2_.push_back(neighbor_node);
            logger_->log(LogLevel::TRACE,"  Adding neighbor {} of initial root to re-rooting stack.", neighbor_node);
        }
    }

    while (!dfs_stack2_.empty()) {
        NodeId current_node = dfs_stack2_.back();
        dfs_stack2_.pop_back();
        logger_->log(LogLevel::TRACE,"  Re-rooting check at node {}", current_node);

        assert(static_cast<size_t>(current_node) < strong_pruning_parent_.size());
        NodeId parent_node = strong_pruning_parent_[current_node].first;
        double parent_edge_cost = strong_pruning_parent_[current_node].second;
        assert(parent_node != kInvalidNodeId);
        assert(static_cast<size_t>(parent_node) < strong_pruning_payoff_.size());
        assert(static_cast<size_t>(current_node) < strong_pruning_payoff_.size());

        double parent_val_without_current = strong_pruning_payoff_[parent_node];
        double current_node_net_payoff = strong_pruning_payoff_[current_node] - parent_edge_cost;

        if (current_node_net_payoff > 1e-9) {
            parent_val_without_current -= current_node_net_payoff;
        }
        logger_->log(LogLevel::TRACE,"    Parent {} payoff: {:.4f}. Current net payoff: {:.4f}. Parent val w/o current: {:.4f}",
                     parent_node, strong_pruning_payoff_[parent_node], current_node_net_payoff, parent_val_without_current);

        if (parent_val_without_current > parent_edge_cost + 1e-9) {
            double contribution_from_parent_side = parent_val_without_current - parent_edge_cost;
            logger_->log(LogLevel::TRACE,"    Adding contribution {:.4f} from parent side to node {}.", contribution_from_parent_side, current_node);
            strong_pruning_payoff_[current_node] += contribution_from_parent_side;
        }

        logger_->log(LogLevel::TRACE,"    Node {} re-rooted payoff: {:.4f}", current_node, strong_pruning_payoff_[current_node]);

        if (strong_pruning_payoff_[current_node] > current_best_value) {
            current_best_root = current_node;
            current_best_value = strong_pruning_payoff_[current_node];
            logger_->log(LogLevel::DEBUG,"    New best root found: {} with payoff {:.4f}", current_best_root, current_best_value);
        }

        assert(static_cast<size_t>(current_node) < neighbors_.size());
        for (const auto& edge_pair : neighbors_[current_node]) {
            NodeId neighbor_node = edge_pair.first;
            if (neighbor_node != parent_node) {

                assert(static_cast<size_t>(neighbor_node) < final_component_label_.size());
                if (final_component_label_[neighbor_node] == component_index) {
                    dfs_stack2_.push_back(neighbor_node);
                    logger_->log(LogLevel::TRACE,"    Adding child {} to re-rooting stack.", neighbor_node);
                }
            }
        }
    }

    logger_->log(LogLevel::TRACE,"Exiting find_best_component_root. Best root: {}", current_best_root);
    return current_best_root;
}

void StrongPruner::mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index) {
    logger_->log(LogLevel::TRACE, "  Marking deleted starting from node {}, parent {}", start_node_index, parent_node_index);
    assert(static_cast<size_t>(start_node_index) < num_nodes_);
    assert(parent_node_index == kInvalidNodeId || static_cast<size_t>(parent_node_index) < num_nodes_);

    node_queue_.clear();

    if (!node_deleted_[start_node_index]) {
        node_deleted_[start_node_index] = true;
        node_queue_.push_back(start_node_index);
        logger_->log(LogLevel::TRACE, "    Marked node {} as deleted.", start_node_index);
    } else {
        logger_->log(LogLevel::TRACE, "    Node {} was already marked deleted.", start_node_index);
        return;
    }

    size_t current_idx = 0;
    while(current_idx < node_queue_.size()) {
        NodeId current_node = node_queue_[current_idx++];
        assert(static_cast<size_t>(current_node) < neighbors_.size());

        for(const auto& edge_pair : neighbors_[current_node]) {
            NodeId neighbor_node = edge_pair.first;

            if (neighbor_node == parent_node_index) {
                continue;
            }

            if (!node_deleted_[neighbor_node]) {
                node_deleted_[neighbor_node] = true;
                node_queue_.push_back(neighbor_node);
                logger_->log(LogLevel::TRACE, "    Marked node {} as deleted (neighbor of {}).", neighbor_node, current_node);
            }
        }

        parent_node_index = kInvalidNodeId;
    }
}

}
}