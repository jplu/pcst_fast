#include "pruning_strategy.h"
// No need to include pcst_fast.h here again, it's in pruning_strategy.h

#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdarg>
#include <memory> // For make_unique
#include <stdexcept> // For invalid_argument in factory

namespace cluster_approx {
namespace internal {

// --- Pruner Implementations ---

// --- NoPruner ---
class NoPruner : public IPruner {
public:
    void prune(const PruningContext& context,
               std::vector<PCSTFast::IndexType>& result_nodes,
               std::vector<PCSTFast::IndexType>& result_edges) override {
        context.log(1, "Pruning: None. Using Phase 1 result directly.\n");
        // Convert phase1_result (vector<int>) to IndexType if necessary,
        // but since IndexType is int, direct copy is fine.
        result_edges.assign(context.phase1_result.begin(), context.phase1_result.end());
        build_node_set(context, result_edges, result_nodes);
        context.log(1, "Final Result (No Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
    }
private:
    void build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
         context.log(3, "NoPruner::build_node_set Entry (using %zu edges).\n", edge_set.size());
         node_set.clear();
         size_t num_nodes = context.prizes.size();
         node_set.reserve(num_nodes);
         std::vector<uint8_t> included_nodes_local(num_nodes, 0);

         for (PCSTFast::IndexType edge_idx : edge_set) { // Use IndexType
             if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;
             const auto& edge = context.edges[edge_idx];
             PCSTFast::IndexType uu = edge.first;
             PCSTFast::IndexType vv = edge.second;
             if (static_cast<size_t>(uu) < included_nodes_local.size() && !included_nodes_local[uu]) {
                 included_nodes_local[uu] = 1; node_set.push_back(uu);
             }
             if (static_cast<size_t>(vv) < included_nodes_local.size() && !included_nodes_local[vv]) {
                 included_nodes_local[vv] = 1; node_set.push_back(vv);
             }
         }
         for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
             bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
             bool is_included = (static_cast<size_t>(ii) < included_nodes_local.size() && included_nodes_local[ii]);
             if (is_good && !is_included) { node_set.push_back(ii); }
         }
         std::sort(node_set.begin(), node_set.end());
         context.log(3, "NoPruner::build_node_set Exit. Final node set size: %zu\n", node_set.size());
    }
};

// --- SimplePruner ---
class SimplePruner : public IPruner {
private:
    std::vector<PCSTFast::IndexType> phase2_result_local_; // Use IndexType

    void build_phase2(const PruningContext& context) {
        context.log(3, "SimplePruner::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
        phase2_result_local_.clear();
        phase2_result_local_.reserve(context.phase1_result.size());
        for (int edge_idx_int : context.phase1_result) { // Iterate original type
            PCSTFast::IndexType edge_idx = static_cast<PCSTFast::IndexType>(edge_idx_int); // Cast if needed
            if (static_cast<size_t>(edge_idx) < context.edges.size()) {
                const auto& edge = context.edges[edge_idx];
                bool u_good = static_cast<size_t>(edge.first) < context.node_good.size() && context.node_good[edge.first];
                bool v_good = static_cast<size_t>(edge.second) < context.node_good.size() && context.node_good[edge.second];
                if (u_good && v_good) {
                    phase2_result_local_.push_back(edge_idx);
                } else {
                     context.log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                }
            }
        }
        context.log(1, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
    }

    void build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
         context.log(3, "SimplePruner::build_node_set Entry (using %zu edges).\n", edge_set.size());
         node_set.clear();
         size_t num_nodes = context.prizes.size();
         node_set.reserve(num_nodes);
         std::vector<uint8_t> included_nodes_local(num_nodes, 0);

         for (PCSTFast::IndexType edge_idx : edge_set) {
              if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;
              const auto& edge = context.edges[edge_idx];
              PCSTFast::IndexType uu = edge.first;
              PCSTFast::IndexType vv = edge.second;
              if (static_cast<size_t>(uu) < included_nodes_local.size() && !included_nodes_local[uu]) {
                  included_nodes_local[uu] = 1; node_set.push_back(uu);
              }
              if (static_cast<size_t>(vv) < included_nodes_local.size() && !included_nodes_local[vv]) {
                  included_nodes_local[vv] = 1; node_set.push_back(vv);
              }
          }
          for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
              bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
              bool is_included = (static_cast<size_t>(ii) < included_nodes_local.size() && included_nodes_local[ii]);
              if (is_good && !is_included) { node_set.push_back(ii); }
          }
          std::sort(node_set.begin(), node_set.end());
          context.log(3, "SimplePruner::build_node_set Exit. Final node set size: %zu\n", node_set.size());
    }

public:
    void prune(const PruningContext& context,
               std::vector<PCSTFast::IndexType>& result_nodes,
               std::vector<PCSTFast::IndexType>& result_edges) override {
        context.log(1, "Pruning: Simple. Running Phase 2 filtering.\n");
        build_phase2(context);
        result_edges = phase2_result_local_; // Copy phase 2 edges
        build_node_set(context, result_edges, result_nodes);
        context.log(1, "Final Result (Simple Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
    }
};


// --- Base class for GW and Strong Pruning (common parts) ---
class AdvancedPrunerBase : public IPruner {
protected:
    // Common state needed by GW/Strong
    std::vector<PCSTFast::IndexType> phase2_result_local_; // Use IndexType
    std::vector<PCSTFast::IndexType> phase3_result_local_; // Use IndexType
    std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>> phase3_neighbors_;
    std::vector<uint8_t> node_deleted_; // Local copy for pruning phase
    std::vector<PCSTFast::IndexType> cluster_queue_; // Local temporary queue

    // Common setup steps
    void build_phase2(const PruningContext& context) {
        context.log(3, "AdvancedPrunerBase::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
        phase2_result_local_.clear();
        phase2_result_local_.reserve(context.phase1_result.size());
        for (int edge_idx_int : context.phase1_result) {
            PCSTFast::IndexType edge_idx = static_cast<PCSTFast::IndexType>(edge_idx_int);
            if (static_cast<size_t>(edge_idx) < context.edges.size()) {
                const auto& edge = context.edges[edge_idx];
                bool u_good = static_cast<size_t>(edge.first) < context.node_good.size() && context.node_good[edge.first];
                bool v_good = static_cast<size_t>(edge.second) < context.node_good.size() && context.node_good[edge.second];
                if (u_good && v_good) {
                    phase2_result_local_.push_back(edge_idx);
                } else {
                     context.log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                }
            }
        }
        context.log(1, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
    }

    void build_phase3_adjacency(const PruningContext& context) {
         context.log(3, "AdvancedPrunerBase::build_phase3_adjacency Entry (Using %zu phase 2 edges).\n", phase2_result_local_.size());
         size_t num_nodes = context.prizes.size();
         if (phase3_neighbors_.size() != num_nodes) phase3_neighbors_.resize(num_nodes);
         for (auto& neighbors : phase3_neighbors_) { neighbors.clear(); }

         for (PCSTFast::IndexType edge_idx : phase2_result_local_) { // Use IndexType
             if (static_cast<size_t>(edge_idx) < context.edges.size() && static_cast<size_t>(edge_idx) < context.costs.size()) {
                 const auto& edge = context.edges[edge_idx];
                 PCSTFast::ValueType cost = context.costs[edge_idx];
                 if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                     phase3_neighbors_[edge.first].emplace_back(edge.second, cost);
                     phase3_neighbors_[edge.second].emplace_back(edge.first, cost);
                 }
             }
         }
          context.log(3, "AdvancedPrunerBase::build_phase3_adjacency Exit.\n");
    }

    // Helper to build final node set based on node_good and local node_deleted_
    void build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set) {
         context.log(3, "AdvancedPrunerBase::build_pruned_node_set Entry.\n");
         node_set.clear();
         size_t num_nodes = context.prizes.size();
         node_set.reserve(num_nodes);
         for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
             bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
             bool is_deleted = (static_cast<size_t>(ii) < node_deleted_.size() && node_deleted_[ii]);
             if (is_good && !is_deleted) {
                 node_set.push_back(ii);
             }
         }
         std::sort(node_set.begin(), node_set.end());
         context.log(3, "AdvancedPrunerBase::build_pruned_node_set Exit. Final node set size: %zu\n", node_set.size());
    }

public:
    virtual ~AdvancedPrunerBase() = default; // Add virtual destructor

    // Common pre-pruning steps for GW/Strong
    virtual void setup(const PruningContext& context) {
        cluster_queue_.reserve(context.prizes.size());
        node_deleted_.assign(context.prizes.size(), 0);
        build_phase2(context);
        build_phase3_adjacency(context);
    }
};


// --- GW Pruner ---
class GWPruner : public AdvancedPrunerBase {
private:
    std::vector<bool> cluster_necessary_local_;

    void mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index) {
        context.log(4, "GWPruner::mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);
         PCSTFast::IndexType current_cluster_index = start_node_index;
         int steps = 0;
         while (current_cluster_index != PCSTFast::kInvalidIndex && static_cast<size_t>(current_cluster_index) < context.clusters.size()) {
             steps++;
             if (static_cast<size_t>(current_cluster_index) < cluster_necessary_local_.size() && cluster_necessary_local_[current_cluster_index]) {
                  context.log(4, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
                 break;
             }
              context.log(4, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);
              if(static_cast<size_t>(current_cluster_index) < cluster_necessary_local_.size()) {
                   cluster_necessary_local_[current_cluster_index] = true;
              } else {
                   context.log(1, "Warning: Cluster index %d out of bounds for necessary flag (%zu).\n", current_cluster_index, cluster_necessary_local_.size());
                   break;
              }
             current_cluster_index = context.clusters[current_cluster_index].merged_into;
              context.log(4, "  Moving up to parent cluster %d.\n", current_cluster_index);
         }
           context.log(4, "GWPruner::mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
    }

    void mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index) {
         // (Implementation unchanged, uses base class members cluster_queue_, node_deleted_)
         context.log(4, "GWPruner::mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);
         if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
              context.log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
             return;
         }
         cluster_queue_.clear();
         cluster_queue_.push_back(start_node_index);
         node_deleted_[start_node_index] = 1;
          context.log(3, "  GW: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);
         size_t queue_index = 0;
         int nodes_deleted_count = 1;
         while (queue_index < cluster_queue_.size()) {
             PCSTFast::IndexType current_node_index = cluster_queue_[queue_index++];
             if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                 for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                     if (neighbor_node_index == parent_node_index) continue;
                     if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                         node_deleted_[neighbor_node_index] = 1;
                         cluster_queue_.push_back(neighbor_node_index);
                         nodes_deleted_count++;
                          context.log(3, "    GW: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                     }
                 }
             }
         }
          context.log(4, "GWPruner::mark_nodes_as_deleted_gw Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
    }

    void run_gw_pruning(const PruningContext& context) {
        context.log(3, "GWPruner::run_gw_pruning Entry.\n");
        phase3_result_local_.clear();
        phase3_result_local_.reserve(phase2_result_local_.size());
        cluster_necessary_local_.assign(context.clusters.size(), false);

        for (int ii = std::ssize(phase2_result_local_) - 1; ii >= 0; --ii) {
            PCSTFast::IndexType edge_idx = phase2_result_local_[ii]; // Use IndexType
             context.log(4, "GW reverse pass: Processing edge %d.\n", edge_idx);

            if(static_cast<size_t>(edge_idx) >= context.edges.size() ||
               static_cast<size_t>(edge_idx) >= context.edge_info.size()) {
                    context.log(1,"Warning: Invalid edge index %d in GW prune loop.\n", edge_idx);
                    continue;
            }

            const auto& edge = context.edges[edge_idx];
            PCSTFast::IndexType uu = edge.first;
            PCSTFast::IndexType vv = edge.second;
            bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
            bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);

            if (u_deleted && v_deleted) {
                 context.log(3, "  GW: Both endpoints (%d, %d) deleted. Skipping edge %d.\n", uu, vv, edge_idx);
                 continue;
            }

            PCSTFast::IndexType inactive_merge_idx = context.edge_info[edge_idx].inactive_merge_event;

            if (inactive_merge_idx == PCSTFast::kInvalidIndex) {
                 context.log(3, "  GW: Edge %d (%d, %d) Active-Active. Keeping.\n", edge_idx, uu, vv);
                 phase3_result_local_.push_back(edge_idx);
                 if (!u_deleted) mark_clusters_as_necessary_gw(context, uu);
                 if (!v_deleted) mark_clusters_as_necessary_gw(context, vv);
            } else {
                 if (static_cast<size_t>(inactive_merge_idx) >= context.inactive_merge_events.size()) {
                      context.log(0,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, edge_idx);
                      continue;
                 }
                 const auto& merge_event = context.inactive_merge_events[inactive_merge_idx];
                 PCSTFast::IndexType inactive_rep = merge_event.inactive_cluster_index;
                 bool is_needed = (static_cast<size_t>(inactive_rep) < cluster_necessary_local_.size() && cluster_necessary_local_[inactive_rep]);

                 if (is_needed) {
                      context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d needed. Keeping.\n", edge_idx, uu, vv, inactive_rep);
                      phase3_result_local_.push_back(edge_idx);
                      mark_clusters_as_necessary_gw(context, merge_event.active_cluster_node);
                      mark_clusters_as_necessary_gw(context, merge_event.inactive_cluster_node);
                 } else {
                      context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d not needed. Pruning from node %d.\n", edge_idx, uu, vv, inactive_rep, merge_event.inactive_cluster_node);
                      mark_nodes_as_deleted_gw(context, merge_event.inactive_cluster_node, merge_event.active_cluster_node);
                 }
            }
        }
        std::reverse(phase3_result_local_.begin(), phase3_result_local_.end());
         context.log(3, "GWPruner::run_gw_pruning Exit. Final edge count: %zu\n", phase3_result_local_.size());
    }

public:
    void prune(const PruningContext& context,
               std::vector<PCSTFast::IndexType>& result_nodes,
               std::vector<PCSTFast::IndexType>& result_edges) override {
        context.log(1, "Pruning: GW. Setting up...\n");
        setup(context);
        context.log(1, "Pruning: Running GW pruning logic...\n");
        run_gw_pruning(context);
        result_edges = phase3_result_local_;
        context.log(2, "GW pruning complete. Building final node set...\n");
        build_pruned_node_set(context, result_nodes); // Use base class helper
        context.log(1, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
    }
};


// --- Strong Pruner ---
class StrongPruner : public AdvancedPrunerBase {
private:
    // Strong pruning specific state
    std::vector<PCSTFast::IndexType> final_component_label_;
    std::vector<std::vector<PCSTFast::IndexType>> final_components_;
    PCSTFast::IndexType root_component_index_ = PCSTFast::kInvalidIndex;
    std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>> strong_pruning_parent_;
    std::vector<PCSTFast::ValueType> strong_pruning_payoff_;
    std::vector<std::pair<bool, PCSTFast::IndexType>> strong_pruning_stack_; // Stack 1
    std::vector<PCSTFast::IndexType> strong_pruning_stack2_; // Stack 2

    // Methods (implementations unchanged from previous fix, just use context.log)
    void label_final_components(const PruningContext& context) {
        context.log(3, "StrongPruner::label_final_components Entry.\n");
         size_t num_nodes = context.prizes.size();
         final_component_label_.assign(num_nodes, PCSTFast::kInvalidIndex);
         final_components_.clear();
         root_component_index_ = PCSTFast::kInvalidIndex;
         int components_found = 0;
         for (PCSTFast::IndexType start_node = 0; start_node < static_cast<PCSTFast::IndexType>(num_nodes); ++start_node) {
              bool is_in_phase2_graph = (static_cast<size_t>(start_node) < phase3_neighbors_.size() && !phase3_neighbors_[start_node].empty())
                                      || (static_cast<size_t>(start_node) < context.node_good.size() && context.node_good[start_node]
                                          && static_cast<size_t>(start_node) < phase3_neighbors_.size() && phase3_neighbors_[start_node].empty());
             if (is_in_phase2_graph && final_component_label_[start_node] == PCSTFast::kInvalidIndex) {
                 PCSTFast::IndexType new_comp_id = final_components_.size();
                 components_found++;
                 context.log(3, "  Found new component %d starting from node %d.\n", new_comp_id, start_node);
                 final_components_.emplace_back();
                 label_component_recursive(context, start_node, new_comp_id);
                 context.log(3, "  Component %d labeled. Size=%zu.\n", new_comp_id, final_components_[new_comp_id].size());
             }
         }
          context.log(3, "StrongPruner::label_final_components Exit. Found %d components.\n", components_found);
    }
    void label_component_recursive(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id) {
         context.log(4, "StrongPruner::label_component_recursive Entry: Start=%d, ID=%d\n", start_node, comp_id);
         cluster_queue_.clear(); // Use base class queue
         cluster_queue_.push_back(start_node);
         final_component_label_[start_node] = comp_id;
         final_components_[comp_id].push_back(start_node);
         if (start_node == context.root) root_component_index_ = comp_id;
         size_t q_idx = 0;
         while(q_idx < cluster_queue_.size()) {
             PCSTFast::IndexType u = cluster_queue_[q_idx++];
             if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                 for(const auto& edge : phase3_neighbors_[u]) {
                     PCSTFast::IndexType v = edge.first;
                     if (static_cast<size_t>(v) < final_component_label_.size() && final_component_label_[v] == PCSTFast::kInvalidIndex) {
                         final_component_label_[v] = comp_id;
                         final_components_[comp_id].push_back(v);
                         cluster_queue_.push_back(v);
                         if (v == context.root) root_component_index_ = comp_id;
                     }
                 }
             }
         }
         context.log(4, "StrongPruner::label_component_recursive Exit for Component %d.\n", comp_id);
    }
     PCSTFast::IndexType find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx) {
          context.log(3, "StrongPruner::find_best_component_root Entry: Component=%d\n", comp_idx);
         if (static_cast<size_t>(comp_idx) >= final_components_.size() || final_components_[comp_idx].empty()) return PCSTFast::kInvalidIndex;
         const auto& comp_nodes = final_components_[comp_idx];
         PCSTFast::IndexType initial_root = comp_nodes[0];
          context.log(4,"  Resetting parent/payoff for component %d nodes.\n", comp_idx);
         for (PCSTFast::IndexType node_idx : comp_nodes) {
             if (static_cast<size_t>(node_idx) < strong_pruning_parent_.size()) strong_pruning_parent_[node_idx] = {PCSTFast::kInvalidIndex, 0.0};
             if (static_cast<size_t>(node_idx) < strong_pruning_payoff_.size()) strong_pruning_payoff_[node_idx] = 0.0;
         }
          context.log(4,"  Running initial strong_pruning_dfs from %d (MarkDeleted=false).\n", initial_root);
         strong_pruning_dfs(context, initial_root, false);
         if (static_cast<size_t>(initial_root) >= strong_pruning_payoff_.size()) return PCSTFast::kInvalidIndex;
         PCSTFast::IndexType best_root = initial_root;
         PCSTFast::ValueType best_value = strong_pruning_payoff_[initial_root];
         context.log(3, "  Initial Payoff at root %d = %.9g\n", best_root, best_value);
          context.log(4,"  Starting payoff propagation from initial root %d.\n", initial_root);
         propagate_payoffs_and_find_best(context, initial_root, best_root, best_value);
          context.log(3, "StrongPruner::find_best_component_root Exit. Best Root=%d, Best Payoff=%.9g\n", best_root, best_value);
         return best_root;
    }
    void propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root, PCSTFast::ValueType& best_value) {
        // (Implementation unchanged, uses context.log and local stack strong_pruning_stack2_)
         context.log(4,"StrongPruner::propagate_payoffs Entry: InitialRoot=%d\n", initial_root);
         strong_pruning_stack2_.clear();
         size_t num_nodes = context.prizes.size();
         if (static_cast<size_t>(initial_root) < phase3_neighbors_.size()) {
              context.log(4, "  Adding children of initial root %d to propagation stack...\n", initial_root);
             for (const auto& edge : phase3_neighbors_[initial_root]) {
                 PCSTFast::IndexType neighbor = edge.first;
                 if (static_cast<size_t>(neighbor) < strong_pruning_parent_.size() && strong_pruning_parent_[neighbor].first == initial_root) {
                      context.log(4, "    Adding child %d to stack.\n", neighbor);
                     strong_pruning_stack2_.push_back(neighbor);
                 }
             }
         }
         int nodes_processed = 0;
         while (!strong_pruning_stack2_.empty()) {
             PCSTFast::IndexType u = strong_pruning_stack2_.back(); strong_pruning_stack2_.pop_back();
             nodes_processed++;
             context.log(4,"    Processing node %d from propagation stack (item %d).\n", u, nodes_processed);
             if (static_cast<size_t>(u) >= num_nodes || static_cast<size_t>(u) >= strong_pruning_parent_.size() || static_cast<size_t>(u) >= strong_pruning_payoff_.size()) continue;
             PCSTFast::IndexType p = strong_pruning_parent_[u].first;
             PCSTFast::ValueType edge_cost = strong_pruning_parent_[u].second;
             context.log(4,"      Parent=%d, EdgeCost=%.9g\n", p, edge_cost);
             if (p == PCSTFast::kInvalidIndex || static_cast<size_t>(p) >= strong_pruning_payoff_.size()) continue;
             PCSTFast::ValueType payoff_u = strong_pruning_payoff_[u];
             PCSTFast::ValueType payoff_p = strong_pruning_payoff_[p];
             context.log(4,"      Payoffs: Subtree(Current)=%.9g, Parent(Rooted)=%.9g\n", payoff_u, payoff_p);
             PCSTFast::ValueType u_contrib = (payoff_u > edge_cost) ? (payoff_u - edge_cost) : 0.0;
             context.log(4,"      Current contribution to parent = max(0, %.9g - %.9g) = %.9g\n", payoff_u, edge_cost, u_contrib);
             PCSTFast::ValueType p_without_u = payoff_p - u_contrib;
             context.log(4,"      Parent payoff without current = %.9g - %.9g = %.9g\n", payoff_p, u_contrib, p_without_u);
             PCSTFast::ValueType p_contrib = (p_without_u > edge_cost) ? (p_without_u - edge_cost) : 0.0;
             context.log(4,"      Parent contribution to current = max(0, %.9g - %.9g) = %.9g\n", p_without_u, edge_cost, p_contrib);
             PCSTFast::ValueType u_total_payoff = payoff_u + p_contrib;
             context.log(4,"      Total payoff if node %d is root = %.9g + %.9g = %.9g\n", u, payoff_u, p_contrib, u_total_payoff);
             if (u_total_payoff > best_value) {
                  context.log(3,"      New best root found: Node=%d, Payoff=%.9g (OldBest: Node=%d, Payoff=%.9g)\n", u, u_total_payoff, best_root, best_value);
                 best_root = u; best_value = u_total_payoff;
             }
             strong_pruning_payoff_[u] = u_total_payoff;
             if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                  context.log(4,"      Adding children of node %d to propagation stack...\n", u);
                 for (const auto& edge : phase3_neighbors_[u]) {
                     PCSTFast::IndexType v = edge.first;
                     if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u) {
                           context.log(4,"        Adding child %d to stack.\n", v);
                         strong_pruning_stack2_.push_back(v);
                     }
                 }
             }
         }
         context.log(4, "StrongPruner::propagate_payoffs Exit. Processed %d nodes.\n", nodes_processed);
    }
    void strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted) {
        // (Implementation unchanged, uses context.log and local stack strong_pruning_stack_)
         context.log(3,"StrongPruner::strong_pruning_dfs Entry: Start=%d, MarkDeleted=%d\n", start_node, mark_deleted);
         strong_pruning_stack_.clear();
         size_t num_nodes = context.prizes.size();
         if(static_cast<size_t>(start_node) >= num_nodes || static_cast<size_t>(start_node) >= strong_pruning_parent_.size()) return;
         context.log(4,"  Initializing DFS root %d: Parent=Invalid, Pushing PreOrder.\n", start_node);
         strong_pruning_parent_[start_node] = {PCSTFast::kInvalidIndex, 0.0};
         strong_pruning_stack_.emplace_back(true, start_node);
         int nodes_visited_pre = 0;
         int nodes_visited_post = 0;
         while(!strong_pruning_stack_.empty()) {
             bool is_pre = strong_pruning_stack_.back().first;
             PCSTFast::IndexType u = strong_pruning_stack_.back().second;
             strong_pruning_stack_.pop_back();
             if(static_cast<size_t>(u) >= num_nodes) continue;
             if(is_pre) {
                 nodes_visited_pre++;
                 context.log(4,"  DFS PreOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_pre);
                 context.log(4,"    Pushing PostOrder visit for node %d.\n", u);
                 strong_pruning_stack_.emplace_back(false, u);
                 strong_pruning_payoff_[u] = context.prizes[u];
                 context.log(4,"    Initialized Payoff[%d] = Prize = %.9g\n", u, strong_pruning_payoff_[u]);
                 if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                      context.log(4,"    Exploring neighbors of %d...\n", u);
                     for(const auto& edge : phase3_neighbors_[u]) {
                         PCSTFast::IndexType v = edge.first;
                         context.log(4,"      Neighbor %d (Parent is %d).\n", v, strong_pruning_parent_[u].first);
                          if(v == strong_pruning_parent_[u].first) { context.log(4,"        Is parent, skipping.\n"); continue; }
                          if(static_cast<size_t>(v) >= num_nodes) { context.log(1,"Warning: Invalid neighbor index %d for node %d.\n", v, u); continue; }
                          context.log(4,"      Setting Parent[%d] = %d (Cost=%.9g). Pushing PreOrder.\n", v, u, edge.second);
                          strong_pruning_parent_[v] = {u, edge.second};
                          strong_pruning_stack_.emplace_back(true, v);
                     }
                 }
             } else {
                 nodes_visited_post++;
                 context.log(4,"  DFS PostOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_post);
                 if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                      context.log(4,"    Aggregating child payoffs for node %d (CurrentPayoff=%.9g)\n", u, strong_pruning_payoff_[u]);
                     for(const auto& edge : phase3_neighbors_[u]) {
                         PCSTFast::IndexType v = edge.first;
                         if(static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u) {
                             context.log(4,"      Processing child %d...\n", v);
                             PCSTFast::ValueType v_cost = strong_pruning_parent_[v].second;
                             PCSTFast::ValueType v_net_payoff = strong_pruning_payoff_[v] - v_cost;
                             context.log(4,"        ChildPayoff=%.9g, EdgeCost=%.9g -> NetChildPayoff=%.9g\n", strong_pruning_payoff_[v], v_cost, v_net_payoff);
                             if(v_net_payoff <= 0.0) {
                                  context.log(3,"        Child %d subtree payoff %.9g <= 0.\n", v, v_net_payoff);
                                 if(mark_deleted) {
                                      context.log(3, "        Strong Pruning: Pruning subtree at %d.\n", v);
                                     mark_nodes_as_deleted_strong(context, v, u);
                                 } else { context.log(4,"        (Not marking deleted as MarkDeleted=false)\n"); }
                             } else {
                                  context.log(4,"        Adding %.9g to Payoff[%d].\n", v_net_payoff, u);
                                 strong_pruning_payoff_[u] += v_net_payoff;
                             }
                         }
                     }
                       context.log(4,"    Final Payoff[%d]=%.9g after aggregating children.\n", u, strong_pruning_payoff_[u]);
                 }
             }
         }
         context.log(3, "StrongPruner::strong_pruning_dfs Exit. Visited %d nodes (pre), %d nodes (post).\n", nodes_visited_pre, nodes_visited_post);
    }
    void mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node) {
        // (Implementation unchanged, uses context.log and base class members)
        context.log(4,"StrongPruner::mark_nodes_as_deleted_strong Entry: Start=%d, Parent=%d\n", start_node, parent_node);
         if (static_cast<size_t>(start_node) >= node_deleted_.size() || node_deleted_[start_node]) { context.log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node); return; }
         cluster_queue_.clear();
         cluster_queue_.push_back(start_node);
         node_deleted_[start_node] = 1;
          context.log(3,"  Strong: Marking node %d and subtree (excluding %d) as deleted.\n", start_node, parent_node);
         size_t q_idx = 0;
         int count = 1;
         while(q_idx < cluster_queue_.size()) {
             PCSTFast::IndexType u = cluster_queue_[q_idx++];
             if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                 for(const auto& edge : phase3_neighbors_[u]) {
                     PCSTFast::IndexType v = edge.first;
                     if(v == parent_node) continue;
                     if(static_cast<size_t>(v) < node_deleted_.size() && !node_deleted_[v]) {
                         node_deleted_[v] = 1;
                         cluster_queue_.push_back(v);
                         count++;
                         context.log(3,"    Strong: Deleted node %d (neighbor of %d).\n", v, u);
                     }
                 }
             }
         }
         context.log(4,"StrongPruner::mark_nodes_as_deleted_strong Exit. Deleted %d nodes.\n", count);
    }

    void run_strong_pruning(const PruningContext& context) {
         context.log(3, "StrongPruner::run_strong_pruning Entry.\n");
         size_t num_nodes = context.prizes.size();
         strong_pruning_parent_.resize(num_nodes);
         strong_pruning_payoff_.resize(num_nodes);
         strong_pruning_stack_.reserve(num_nodes);
         strong_pruning_stack2_.reserve(num_nodes);

         label_final_components(context);
         context.log(2, "Strong Pruning: Found %zu components.\n", final_components_.size());

         for (size_t comp_idx = 0; comp_idx < final_components_.size(); ++comp_idx) {
             if(final_components_[comp_idx].empty()) continue;
             context.log(2,"Strong Pruning: Processing component %zu (size %zu).\n", comp_idx, final_components_[comp_idx].size());
             PCSTFast::IndexType root_node;
             if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                 root_node = context.root;
                 context.log(3,"  Using designated root %d.\n", root_node);
             } else {
                 root_node = find_best_component_root(context, comp_idx);
                 context.log(3,"  Using best root %d.\n", root_node);
             }

             if(root_node != PCSTFast::kInvalidIndex && static_cast<size_t>(root_node) < num_nodes) {
                 for(PCSTFast::IndexType n : final_components_[comp_idx]) {
                     if(static_cast<size_t>(n) < strong_pruning_parent_.size()) strong_pruning_parent_[n] = {PCSTFast::kInvalidIndex, 0.0};
                     if(static_cast<size_t>(n) < strong_pruning_payoff_.size()) strong_pruning_payoff_[n] = 0.0;
                 }
                 strong_pruning_dfs(context, root_node, true);
             } else {
                  context.log(1,"Warning: Skipping strong pruning for component %zu due to invalid root.\n", comp_idx);
             }
         }
         // Filter phase 2 edges based on local node_deleted_
         phase3_result_local_.clear();
         phase3_result_local_.reserve(phase2_result_local_.size());
         context.log(2, "Strong pruning complete. Filtering edges based on deleted nodes...\n");
         for(PCSTFast::IndexType edge_idx : phase2_result_local_) { // Use IndexType
             if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;
             const auto& edge = context.edges[edge_idx];
             bool u_del = (static_cast<size_t>(edge.first) >= node_deleted_.size() || node_deleted_[edge.first]);
             bool v_del = (static_cast<size_t>(edge.second) >= node_deleted_.size() || node_deleted_[edge.second]);
             if (!u_del && !v_del) {
                 phase3_result_local_.push_back(edge_idx);
             } else {
                  context.log(4,"  Strong: Edge %d (%d,%d) removed due to deleted endpoint.\n", edge_idx, edge.first, edge.second);
             }
         }
         context.log(3, "StrongPruner::run_strong_pruning Exit. Final edge count %zu.\n", phase3_result_local_.size());
    }


public:
    void prune(const PruningContext& context,
               std::vector<PCSTFast::IndexType>& result_nodes,
               std::vector<PCSTFast::IndexType>& result_edges) override {
        context.log(1, "Pruning: Strong. Setting up...\n");
        setup(context);
        context.log(1, "Pruning: Running Strong pruning logic...\n");
        run_strong_pruning(context);
        result_edges = phase3_result_local_;
        context.log(2, "Strong pruning complete. Building final node set...\n");
        build_pruned_node_set(context, result_nodes); // Use base class helper
        context.log(1, "Final Result (Strong Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
    }

};


// --- Factory Function Implementation ---
std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method) {
    switch (method) {
        case PCSTFast::PruningMethod::kNoPruning:
            return std::make_unique<NoPruner>();
        case PCSTFast::PruningMethod::kSimplePruning:
            return std::make_unique<SimplePruner>();
        case PCSTFast::PruningMethod::kGWPruning:
            return std::make_unique<GWPruner>();
        case PCSTFast::PruningMethod::kStrongPruning:
            return std::make_unique<StrongPruner>();
        default: // kUnknownPruning or others
            throw std::invalid_argument("Unsupported pruning method provided to factory.");
    }
}


} // namespace internal
} // namespace cluster_approx
