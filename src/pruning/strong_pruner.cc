#include "pcst_fast/pruning/strong_pruner.h"
#include "pcst_fast/pruning/pruning_utils.h"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

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
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;

        if (input.core_result.initial_node_filter[u] && input.core_result.initial_node_filter[v]) {
            intermediate_edges.push_back(edge_idx);
        }
    }

    if (intermediate_edges.empty()) {
        std::vector<uint8_t> temp_deleted(num_nodes_, 0);
        return { build_final_node_set(num_nodes_, temp_deleted, input.core_result.initial_node_filter), {} };
    }

    node_deleted_.assign(num_nodes_, 0);
    final_component_label_.assign(num_nodes_, kInvalidClusterId);
    final_components_.clear();
    root_component_index_ = kInvalidClusterId;
    
    strong_pruning_parent_.assign(num_nodes_, {kInvalidNodeId, 0.0});
    strong_pruning_payoff_.assign(num_nodes_, -1.0);
    
    neighbors_ = build_adjacency_list_csr(intermediate_edges, input.graph);
    std::vector<NodeId> global_dfs_stack2;

    for (NodeId i = 0; i < static_cast<NodeId>(num_nodes_); ++i) {
        uint32_t csr_start = neighbors_.row_ptr[i];
        uint32_t csr_end = neighbors_.row_ptr[i + 1];
        bool has_edges = (csr_start != csr_end);

        if ((has_edges || (static_cast<size_t>(i) < input.core_result.initial_node_filter.size() && input.core_result.initial_node_filter[i])) && 
            final_component_label_[i] == kInvalidClusterId) {
            final_components_.emplace_back();
            ClusterId current_component_idx = final_components_.size() - 1;
            label_final_component(i, current_component_idx, global_dfs_stack2);
        }
    }

    std::vector<ClusterId> valid_components;
    valid_components.reserve(final_components_.size());
    for (ClusterId comp_idx = 0; comp_idx < static_cast<ClusterId>(final_components_.size()); ++comp_idx) {
        if (!final_components_[comp_idx].empty()) valid_components.push_back(comp_idx);
    }

    // Process components fully in parallel. Node ownership is strictly disjoint per thread loop structure avoiding races entirely.
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < static_cast<int>(valid_components.size()); ++i) {
        ClusterId comp_idx = valid_components[i];
        
        std::vector<std::pair<bool, NodeId>> local_dfs_stack;
        std::vector<NodeId> local_dfs_stack2;
        std::vector<NodeId> local_node_queue;

        if (comp_idx == root_component_index_) {
            assert(input_->graph.root != kInvalidNodeId);
            for (NodeId node : final_components_[comp_idx]) {
                strong_pruning_parent_[node] = {kInvalidNodeId, 0.0};
                strong_pruning_payoff_[node] = -1.0;
            }
            strong_pruning_dfs(input_->graph.root, true, local_dfs_stack, local_node_queue);
        } else {
            NodeId best_root = find_best_component_root(comp_idx, local_dfs_stack, local_dfs_stack2);
            for (NodeId node : final_components_[comp_idx]) {
                strong_pruning_parent_[node] = {kInvalidNodeId, 0.0};
                strong_pruning_payoff_[node] = -1.0;
            }
            strong_pruning_dfs(best_root, true, local_dfs_stack, local_node_queue);
        }
    }

    std::vector<EdgeId> final_edges;
    final_edges.reserve(intermediate_edges.size());
    for (EdgeId edge_idx : intermediate_edges) {
        NodeId u = input.graph.edges[edge_idx].first;
        NodeId v = input.graph.edges[edge_idx].second;
        if (!node_deleted_[u] && !node_deleted_[v]) {
            final_edges.push_back(edge_idx);
        }
    }

    PruningResult result;
    result.edges = std::move(final_edges);
    result.nodes = build_final_node_set(num_nodes_, node_deleted_, input.core_result.initial_node_filter);

    input_ = nullptr;
    logger_ = nullptr;

    return result;
}

void StrongPruner::label_final_component(NodeId start_node_index, ClusterId component_index, std::vector<NodeId>& dfs_stack2) {
    dfs_stack2.clear();
    dfs_stack2.push_back(start_node_index);
    final_component_label_[start_node_index] = component_index;

    while (!dfs_stack2.empty()) {
        NodeId current_node = dfs_stack2.back();
        dfs_stack2.pop_back();

        final_components_[component_index].push_back(current_node);
        if (current_node == input_->graph.root) {
            root_component_index_ = component_index;
        }

        for (const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
            NodeId neighbor_node = edge_pair.first;
            if (final_component_label_[neighbor_node] == kInvalidClusterId) {
                final_component_label_[neighbor_node] = component_index;
                dfs_stack2.push_back(neighbor_node);
            }
        }
    }
}

void StrongPruner::strong_pruning_dfs(NodeId start_node_index, bool mark_as_deleted,
                                      std::vector<std::pair<bool, NodeId>>& local_dfs_stack,
                                      std::vector<NodeId>& local_node_queue) {
    local_dfs_stack.clear();
    strong_pruning_parent_[start_node_index] = {kInvalidNodeId, 0.0};
    local_dfs_stack.push_back({true, start_node_index});

    while (!local_dfs_stack.empty()) {
        auto [is_entry_call, current_node] = local_dfs_stack.back();
        local_dfs_stack.pop_back();

        if (is_entry_call) {
            local_dfs_stack.push_back({false, current_node});
            for (const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
                NodeId neighbor_node = edge_pair.first;
                double edge_cost = edge_pair.second;

                if (neighbor_node == strong_pruning_parent_[current_node].first) continue;
                strong_pruning_parent_[neighbor_node] = {current_node, edge_cost};
                local_dfs_stack.push_back({true, neighbor_node});
            }
        } else {
            strong_pruning_payoff_[current_node] = input_->graph.prizes[current_node];

            for (const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
                NodeId neighbor_node = edge_pair.first;
                double edge_cost = edge_pair.second;

                if (strong_pruning_parent_[neighbor_node].first != current_node) continue;

                double child_net_payoff = strong_pruning_payoff_[neighbor_node] - edge_cost;

                if (child_net_payoff <= 1e-9) {
                    if (mark_as_deleted) {
                        mark_nodes_as_deleted(neighbor_node, current_node, local_node_queue);
                    }
                } else {
                    strong_pruning_payoff_[current_node] += child_net_payoff;
                }
            }
        }
    }
}

NodeId StrongPruner::find_best_component_root(ClusterId component_index,
                                              std::vector<std::pair<bool, NodeId>>& local_dfs_stack,
                                              std::vector<NodeId>& local_dfs_stack2) {
    const auto& component_nodes = final_components_[component_index];
    NodeId initial_root = component_nodes[0];

    for (NodeId node : component_nodes) {
        strong_pruning_parent_[node] = {kInvalidNodeId, 0.0};
        strong_pruning_payoff_[node] = -1.0;
    }
    
    std::vector<NodeId> dummy_node_queue;
    strong_pruning_dfs(initial_root, false, local_dfs_stack, dummy_node_queue);

    NodeId current_best_root = initial_root;
    double current_best_value = strong_pruning_payoff_[initial_root];

    local_dfs_stack2.clear();

    for (const auto& edge_pair : neighbors_.get_neighbors(initial_root)) {
        NodeId neighbor_node = edge_pair.first;
        if (final_component_label_[neighbor_node] == component_index) {
            local_dfs_stack2.push_back(neighbor_node);
        }
    }

    while (!local_dfs_stack2.empty()) {
        NodeId current_node = local_dfs_stack2.back();
        local_dfs_stack2.pop_back();

        NodeId parent_node = strong_pruning_parent_[current_node].first;
        double parent_edge_cost = strong_pruning_parent_[current_node].second;

        double parent_val_without_current = strong_pruning_payoff_[parent_node];
        double current_node_net_payoff = strong_pruning_payoff_[current_node] - parent_edge_cost;

        if (current_node_net_payoff > 1e-9) {
            parent_val_without_current -= current_node_net_payoff;
        }

        if (parent_val_without_current > parent_edge_cost + 1e-9) {
            double contribution_from_parent_side = parent_val_without_current - parent_edge_cost;
            strong_pruning_payoff_[current_node] += contribution_from_parent_side;
        }

        if (strong_pruning_payoff_[current_node] > current_best_value) {
            current_best_root = current_node;
            current_best_value = strong_pruning_payoff_[current_node];
        }

        for (const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
            NodeId neighbor_node = edge_pair.first;
            if (neighbor_node != parent_node && final_component_label_[neighbor_node] == component_index) {
                local_dfs_stack2.push_back(neighbor_node);
            }
        }
    }

    return current_best_root;
}

void StrongPruner::mark_nodes_as_deleted(NodeId start_node_index, NodeId parent_node_index,
                                         std::vector<NodeId>& local_node_queue) {
    local_node_queue.clear();

    if (!node_deleted_[start_node_index]) {
        node_deleted_[start_node_index] = 1;
        local_node_queue.push_back(start_node_index);
    } else {
        return;
    }

    size_t current_idx = 0;
    while(current_idx < local_node_queue.size()) {
        NodeId current_node = local_node_queue[current_idx++];

        for(const auto& edge_pair : neighbors_.get_neighbors(current_node)) {
            NodeId neighbor_node = edge_pair.first;

            if (neighbor_node == parent_node_index) continue;

            if (!node_deleted_[neighbor_node]) {
                node_deleted_[neighbor_node] = 1;
                local_node_queue.push_back(neighbor_node);
            }
        }
        parent_node_index = kInvalidNodeId;
    }
}

}
}