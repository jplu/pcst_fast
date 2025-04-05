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
#include <cassert> // For assertions

namespace cluster_approx {
namespace internal {

// --- Pruner Implementations ---

// --- NoPruner ---
void NoPruner::prune(const PruningContext& context,
               std::vector<PCSTFast::IndexType>& result_nodes,
               std::vector<PCSTFast::IndexType>& result_edges) {
    context.log(1, "Pruning: None. Using Phase 1 result directly.\n");
    // Phase1_result contains int indices, result_edges expects IndexType (which is int). Direct assign is okay.
    result_edges.assign(context.phase1_result.begin(), context.phase1_result.end());
    // Build the corresponding node set
    build_node_set(context, result_edges, result_nodes);
    context.log(1, "Final Result (No Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
}

// Helper to build node set from edges and isolated good nodes
void NoPruner::build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
     context.log(3, "NoPruner::build_node_set Entry (using %zu edges).\n", edge_set.size());
     node_set.clear();
     size_t num_nodes = context.prizes.size();
     node_set.reserve(num_nodes);
     // Use a local buffer to track nodes included via edges
     std::vector<uint8_t> included_nodes_local(num_nodes, 0);

     // Add nodes connected by the edges
     for (PCSTFast::IndexType edge_idx : edge_set) { // Use IndexType
         // Basic bounds check
         if (static_cast<size_t>(edge_idx) >= context.edges.size()) {
              context.log(1, "Warning: Invalid edge index %d in NoPruner::build_node_set.\n", edge_idx);
              continue;
         }
         const auto& edge = context.edges[edge_idx];
         PCSTFast::IndexType uu = edge.first;
         PCSTFast::IndexType vv = edge.second;
         // Check bounds before accessing local buffer
         if (static_cast<size_t>(uu) < included_nodes_local.size() && !included_nodes_local[uu]) {
             included_nodes_local[uu] = 1; node_set.push_back(uu);
         }
         if (static_cast<size_t>(vv) < included_nodes_local.size() && !included_nodes_local[vv]) {
             included_nodes_local[vv] = 1; node_set.push_back(vv);
         }
     }

     // Add any 'good' nodes that were isolated (not covered by edges)
     for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
         bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
         bool is_included = (static_cast<size_t>(ii) < included_nodes_local.size() && included_nodes_local[ii]);
         if (is_good && !is_included) {
              context.log(4, "  Adding isolated good node %d.\n", ii);
              node_set.push_back(ii);
         }
     }
     // Sort the final node list
     std::sort(node_set.begin(), node_set.end());
     context.log(3, "NoPruner::build_node_set Exit. Final node set size: %zu\n", node_set.size());
}


// --- SimplePruner ---
void SimplePruner::build_phase2(const PruningContext& context) {
    context.log(3, "SimplePruner::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
    phase2_result_local_.clear(); // Clear previous result if any
    phase2_result_local_.reserve(context.phase1_result.size());
    // Iterate through phase 1 edges (which have int indices)
    for (int edge_idx_int : context.phase1_result) {
        PCSTFast::IndexType edge_idx = static_cast<PCSTFast::IndexType>(edge_idx_int); // Cast to IndexType
        // Bounds check edge index
        if (static_cast<size_t>(edge_idx) < context.edges.size()) {
            const auto& edge = context.edges[edge_idx];
            // Check if both endpoints are 'good' nodes
            bool u_good = static_cast<size_t>(edge.first) < context.node_good.size() && context.node_good[edge.first];
            bool v_good = static_cast<size_t>(edge.second) < context.node_good.size() && context.node_good[edge.second];
            if (u_good && v_good) {
                // Keep edge if both endpoints are connected to the main structure
                phase2_result_local_.push_back(edge_idx);
            } else {
                 context.log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
            }
        } else {
            context.log(1, "Warning: Invalid edge index %d in SimplePruner::build_phase2.\n", edge_idx);
        }
    }
    context.log(1, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
}

// SimplePruner uses the same node building logic as NoPruner for its phase 2 result
void SimplePruner::build_node_set(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
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

void SimplePruner::prune(const PruningContext& context,
           std::vector<PCSTFast::IndexType>& result_nodes,
           std::vector<PCSTFast::IndexType>& result_edges) {
    context.log(1, "Pruning: Simple. Running Phase 2 filtering.\n");
    // Calculate the phase 2 edge set
    build_phase2(context);
    // Assign the phase 2 edges as the final result
    result_edges = phase2_result_local_; // Copy phase 2 edges
    // Build the corresponding node set
    build_node_set(context, result_edges, result_nodes);
    context.log(1, "Final Result (Simple Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
}


// --- Base class for GW and Strong Pruning (common parts) ---
void AdvancedPrunerBase::build_phase2(const PruningContext& context) {
    // Identical logic to SimplePruner::build_phase2
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
        } else {
             context.log(1, "Warning: Invalid edge index %d in AdvancedPrunerBase::build_phase2.\n", edge_idx);
        }
    }
    context.log(1, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
}

void AdvancedPrunerBase::build_phase3_adjacency(const PruningContext& context) {
     context.log(3, "AdvancedPrunerBase::build_phase3_adjacency Entry (Using %zu phase 2 edges).\n", phase2_result_local_.size());
     size_t num_nodes = context.prizes.size();
     // Ensure adjacency list size matches number of nodes
     if (phase3_neighbors_.size() != num_nodes) {
        phase3_neighbors_.resize(num_nodes);
     }
     // Clear lists from previous runs (if any, though typically one run per object)
     for (auto& neighbors : phase3_neighbors_) { neighbors.clear(); }

     int edges_added = 0;
     // Build adjacency list from phase 2 edges
     for (PCSTFast::IndexType edge_idx : phase2_result_local_) {
         // Check edge index validity for accessing edges and costs
         if (static_cast<size_t>(edge_idx) < context.edges.size() && static_cast<size_t>(edge_idx) < context.costs.size()) {
             const auto& edge = context.edges[edge_idx];
             PCSTFast::ValueType cost = context.costs[edge_idx];
             // Check node indices validity before adding to adjacency list
             if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                 // Add edge in both directions
                 phase3_neighbors_[edge.first].emplace_back(edge.second, cost);
                 phase3_neighbors_[edge.second].emplace_back(edge.first, cost);
                 edges_added++;
             } else {
                  context.log(0, "Error: Invalid node index in edge %d while building phase 3 adjacency list.\n", edge_idx);
                  assert(false && "Invalid node index in build_phase3_adjacency");
             }
         } else {
             context.log(1, "Warning: Invalid edge index %d found in phase2_result_ during adjacency build.\n", edge_idx);
         }
     }
      context.log(3, "AdvancedPrunerBase::build_phase3_adjacency Exit. Added %d edges (x2) to lists.\n", edges_added);
}

// Builds the final node set based on 'good' nodes and nodes deleted by the specific pruning strategy
void AdvancedPrunerBase::build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set) {
     context.log(3, "AdvancedPrunerBase::build_pruned_node_set Entry.\n");
     node_set.clear();
     size_t num_nodes = context.prizes.size();
     node_set.reserve(num_nodes);
     int nodes_included = 0;
     // Iterate through all potential nodes
     for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
         // Check if node was identified as 'good' initially
         bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
         // Check if node was deleted by the specific pruning phase (using local node_deleted_)
         bool is_deleted = (static_cast<size_t>(ii) < node_deleted_.size() && node_deleted_[ii]);
         // Include node if it was good AND was not deleted
         if (is_good && !is_deleted) {
             node_set.push_back(ii);
             nodes_included++;
         }
     }
     // Sort the final node list
     std::sort(node_set.begin(), node_set.end());
     context.log(3, "AdvancedPrunerBase::build_pruned_node_set Exit. Final node set size: %zu\n", nodes_included);
}

void AdvancedPrunerBase::setup(const PruningContext& context) {
    context.log(4, "AdvancedPrunerBase::setup Entry.\n");
    // Reserve space for members that will be populated
    phase2_result_local_.reserve(context.phase1_result.size());
    // Resize members based on input size
    size_t num_nodes = context.prizes.size();
    node_deleted_.assign(num_nodes, 0);
    phase3_neighbors_.resize(num_nodes); // Resized again just in case

    // Run common build steps
    build_phase2(context);
    build_phase3_adjacency(context);
    context.log(4, "AdvancedPrunerBase::setup Exit.\n");
}


// --- GW Pruner ---
void GWPruner::mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index) {
    context.log(4, "GWPruner::mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);
    PCSTFast::IndexType current_cluster_index = start_node_index; // Start at the original node index
    int steps = 0;
    // Traverse up the merge tree using merged_into pointers
    while (current_cluster_index != PCSTFast::kInvalidIndex && static_cast<size_t>(current_cluster_index) < context.clusters.size()) {
        steps++;
        // Check bounds for local flag vector
        if (static_cast<size_t>(current_cluster_index) >= cluster_necessary_local_.size()) {
             context.log(1, "Warning: Cluster index %d out of bounds for necessary flag (%zu) during traversal.\n", current_cluster_index, cluster_necessary_local_.size());
             break; // Stop if index is invalid
        }
        // Stop if this cluster (or an ancestor) is already marked
        if (cluster_necessary_local_[current_cluster_index]) {
            context.log(4, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
            break;
        }
        // Mark the current cluster as necessary
        context.log(4, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);
        cluster_necessary_local_[current_cluster_index] = true;
        // Move to the parent cluster
        current_cluster_index = context.clusters[current_cluster_index].merged_into;
        context.log(4, "  Moving up to parent cluster %d.\n", current_cluster_index);
    }
    context.log(4, "GWPruner::mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
}

void GWPruner::mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index) {
    context.log(4, "GWPruner::mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);
    // Check if start node is already deleted or invalid
    if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
        context.log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
        return;
    }

    // Use a local queue for the BFS/DFS traversal
    std::vector<PCSTFast::IndexType> cluster_queue_local;
    cluster_queue_local.reserve(context.prizes.size()); // Estimate capacity

    cluster_queue_local.push_back(start_node_index); // Start traversal from the given node
    node_deleted_[start_node_index] = 1; // Mark the start node as deleted (modifies member)
    context.log(3, "  GW: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);

    size_t queue_index = 0;
    int nodes_deleted_count = 1;
    // Process the queue
    while (queue_index < cluster_queue_local.size()) {
        PCSTFast::IndexType current_node_index = cluster_queue_local[queue_index++];
        context.log(4,"    Processing node %d from GW deletion queue (index %zu).\n", current_node_index, queue_index-1);

        // Explore neighbors in the phase 3 graph (built from phase 2 edges)
        if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
            context.log(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
            for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                context.log(4,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);
                // Don't traverse back towards the 'parent' side which should be kept
                if (neighbor_node_index == parent_node_index) {
                    context.log(4,"          Neighbor is parent, skipping.\n");
                    continue;
                }
                // Check if neighbor is valid and not already deleted
                if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                    node_deleted_[neighbor_node_index] = 1; // Mark neighbor deleted (modify member)
                    cluster_queue_local.push_back(neighbor_node_index); // Add to local queue
                    nodes_deleted_count++;
                    context.log(3, "    GW: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                } else {
                     context.log(4,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
                }
            }
        } else {
             context.log(4, "      Node %d has no neighbors in phase 3 graph (or index out of bounds).\n", current_node_index);
        }
    }
    context.log(4, "GWPruner::mark_nodes_as_deleted_gw Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
}

void GWPruner::run_gw_pruning(const PruningContext& context) {
    context.log(3, "GWPruner::run_gw_pruning Entry.\n");
    phase3_result_local_.clear(); // Clear final edge set
    phase3_result_local_.reserve(phase2_result_local_.size()); // Reserve capacity
    // Resize and reset cluster necessary flags for this run
    cluster_necessary_local_.assign(context.clusters.size(), false);

    context.log(2, "Starting GW pruning reverse pass (processing %zu edges).\n", phase2_result_local_.size());
    int edges_kept = 0;
    int edges_discarded = 0;

    // Process phase 2 edges in reverse order of addition
    for (int ii = std::ssize(phase2_result_local_) - 1; ii >= 0; --ii) {
        PCSTFast::IndexType edge_idx = phase2_result_local_[ii]; // Use IndexType
         context.log(4, "GW reverse pass: Processing edge %d (index %d from end).\n", edge_idx, (int)phase2_result_local_.size()-1-ii);

        // Validate edge index against context data
        if(static_cast<size_t>(edge_idx) >= context.edges.size() ||
           static_cast<size_t>(edge_idx) >= context.edge_info.size()) {
                context.log(1,"Warning: Invalid edge index %d in GW prune loop.\n", edge_idx);
                continue;
        }

        const auto& edge = context.edges[edge_idx];
        PCSTFast::IndexType uu = edge.first;
        PCSTFast::IndexType vv = edge.second;
         context.log(4, "  Edge %d connects nodes (%d, %d).\n", edge_idx, uu, vv);

        // Check if endpoints have already been deleted by *this* pruning phase
        bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
        bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);
         context.log(4, "  Node status: %d Deleted=%d, %d Deleted=%d.\n", uu, u_deleted, vv, v_deleted);

        // If both endpoints are deleted, the edge is implicitly removed
        if (u_deleted && v_deleted) {
             context.log(3, "  GW: Both endpoints (%d, %d) of edge %d deleted. Skipping.\n", uu, vv, edge_idx);
             edges_discarded++;
             continue;
        }

        // Check if this edge caused an Active-Inactive merge
        PCSTFast::IndexType inactive_merge_idx = context.edge_info[edge_idx].inactive_merge_event;
         context.log(4, "  Inactive merge event index for edge %d: %d.\n", edge_idx, inactive_merge_idx);


        if (inactive_merge_idx == PCSTFast::kInvalidIndex) {
             // --- Case 1: Edge from an Active-Active merge ---
             context.log(3, "  GW: Edge %d (%d, %d) Active-Active. Keeping.\n", edge_idx, uu, vv);
             phase3_result_local_.push_back(edge_idx); // Keep the edge
             edges_kept++;
             // Mark the original node clusters as necessary if endpoints not deleted
              context.log(4, "  Marking clusters necessary from nodes %d and %d (if not deleted).\n", uu, vv);
             if (!u_deleted) mark_clusters_as_necessary_gw(context, uu);
             if (!v_deleted) mark_clusters_as_necessary_gw(context, vv);
        } else {
             // --- Case 2: Edge from an Active-Inactive merge ---
             // Validate the inactive merge event index
             if (static_cast<size_t>(inactive_merge_idx) >= context.inactive_merge_events.size()) {
                  context.log(0,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, edge_idx);
                  assert(false && "Invalid inactive merge index");
                  continue;
             }
             const auto& merge_event = context.inactive_merge_events[inactive_merge_idx];
             PCSTFast::IndexType inactive_rep = merge_event.inactive_cluster_index; // Representative of inactive side at merge time
             PCSTFast::IndexType active_node = merge_event.active_cluster_node;    // Original node on active side
             PCSTFast::IndexType inactive_node = merge_event.inactive_cluster_node; // Original node on inactive side
              context.log(4, "  Edge %d was Active-Inactive. ActiveNode=%d, InactiveNode=%d, InactiveRep=%d.\n",
                          edge_idx, active_node, inactive_node, inactive_rep);

             // Check if the inactive representative cluster became necessary later (by being on a path from another kept edge)
             bool is_inactive_necessary = (static_cast<size_t>(inactive_rep) < cluster_necessary_local_.size() && cluster_necessary_local_[inactive_rep]);
              context.log(4, "  Checking necessity of inactive representative cluster %d: Necessary=%d\n",
                          inactive_rep, is_inactive_necessary);

             if (is_inactive_necessary) {
                  // If inactive side is necessary, keep this edge and mark both sides necessary now
                  context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d needed. Keeping.\n", edge_idx, uu, vv, inactive_rep);
                  phase3_result_local_.push_back(edge_idx);
                  edges_kept++;
                   context.log(4, "  Marking clusters necessary from nodes %d and %d.\n", active_node, inactive_node);
                  mark_clusters_as_necessary_gw(context, active_node);
                  mark_clusters_as_necessary_gw(context, inactive_node);
             } else {
                  // If inactive side is not necessary, prune the subtree rooted at the inactive node
                  context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d not needed. Pruning from node %d.\n", edge_idx, uu, vv, inactive_rep, merge_event.inactive_cluster_node);
                  edges_discarded++;
                  // Mark the inactive subtree as deleted, avoiding path back towards active side
                  mark_nodes_as_deleted_gw(context, inactive_node, active_node);
             }
        }
    }

    // Reverse the result as edges were added in reverse order
    context.log(2, "GW pruning reverse pass complete. Reversing phase 3 result.\n");
    std::reverse(phase3_result_local_.begin(), phase3_result_local_.end());
    context.log(3, "GWPruner::run_gw_pruning Exit. Final edge count: %zu (Kept: %d, Discarded: %d)\n", phase3_result_local_.size(), edges_kept, edges_discarded);
}

void GWPruner::prune(const PruningContext& context,
           std::vector<PCSTFast::IndexType>& result_nodes,
           std::vector<PCSTFast::IndexType>& result_edges) {
    context.log(1, "Pruning: GW. Setting up...\n");
    setup(context); // Call base setup (builds phase 2 result, adjacency list)
    // Resize local flag vector based on actual cluster count from context
    cluster_necessary_local_.resize(context.clusters.size());
    context.log(1, "Pruning: Running GW pruning logic...\n");
    run_gw_pruning(context); // Run the specific GW logic
    result_edges = phase3_result_local_; // Assign final GW edges
    context.log(2, "GW pruning complete. Building final node set...\n");
    build_pruned_node_set(context, result_nodes); // Build nodes based on final edges and deletion flags
    context.log(1, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
}


// --- Strong Pruner ---
void StrongPruner::label_final_components(const PruningContext& context) {
    context.log(3, "StrongPruner::label_final_components Entry.\n");
    size_t num_nodes = context.prizes.size();
    // Initialize component tracking members
    final_component_label_.assign(num_nodes, PCSTFast::kInvalidIndex);
    final_components_.clear();
    root_component_index_ = PCSTFast::kInvalidIndex; // Reset root component tracker
    int components_found = 0;

    // Iterate through all nodes to find starting points for unlabeled components
    for (PCSTFast::IndexType start_node = 0; start_node < static_cast<PCSTFast::IndexType>(num_nodes); ++start_node) {
        // Check if this node should be considered (part of phase 2 graph and not yet labeled)
        bool is_in_phase2_graph = (static_cast<size_t>(start_node) < phase3_neighbors_.size() && !phase3_neighbors_[start_node].empty()) // Has neighbors
                                || (static_cast<size_t>(start_node) < context.node_good.size() && context.node_good[start_node]          // Is isolated but good
                                    && static_cast<size_t>(start_node) < phase3_neighbors_.size() && phase3_neighbors_[start_node].empty());

        if (is_in_phase2_graph && final_component_label_[start_node] == PCSTFast::kInvalidIndex) {
            // Found the start of a new component
            PCSTFast::IndexType new_component_id = static_cast<PCSTFast::IndexType>(final_components_.size());
            components_found++;
            context.log(3, "  Found new component %d starting from node %d.\n", new_component_id, start_node);
            final_components_.emplace_back(); // Add space for the new component's node list
            // Recursively (or iteratively) label all nodes connected to start_node
            label_component_recursive(context, start_node, new_component_id);
            context.log(3, "  Component %d labeled. Size=%zu.\n", new_component_id, final_components_[new_component_id].size());
        }
    }
    context.log(3, "StrongPruner::label_final_components Exit. Found %d components.\n", components_found);
}

void StrongPruner::label_component_recursive(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id) {
    context.log(4, "StrongPruner::label_component_recursive Entry: Start=%d, ID=%d\n", start_node, comp_id);
    // Use a local queue for the BFS traversal
    std::vector<PCSTFast::IndexType> cluster_queue_local;
    cluster_queue_local.reserve(context.prizes.size()); // Reserve capacity

    cluster_queue_local.push_back(start_node);
    // Check bounds before accessing members
    assert(static_cast<size_t>(start_node) < final_component_label_.size() && "Start node index out of bounds for label");
    assert(static_cast<size_t>(comp_id) < final_components_.size() && "Component ID out of bounds for components vector");
    final_component_label_[start_node] = comp_id; // Label start node
    final_components_[comp_id].push_back(start_node); // Add start node to component list
    context.log(4, "  Added start node %d to component %d.\n", start_node, comp_id);
    // If this node is the designated root, record its component ID
    if (start_node == context.root) {
        root_component_index_ = comp_id;
        context.log(4, "  Start node is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
    }

    size_t q_idx = 0;
    // Perform BFS
    while(q_idx < cluster_queue_local.size()) {
        PCSTFast::IndexType u = cluster_queue_local[q_idx++];
        context.log(4,"    Processing node %d from component labeling queue (index %zu).\n", u, q_idx-1);

        // Explore neighbors using the phase 3 adjacency list
        if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
            context.log(4,"      Neighbors of %d: %zu\n", u, phase3_neighbors_[u].size());
            for(const auto& edge : phase3_neighbors_[u]) {
                PCSTFast::IndexType v = edge.first; // Neighbor node index
                // Check if neighbor is valid and unlabeled
                if (static_cast<size_t>(v) < final_component_label_.size() && final_component_label_[v] == PCSTFast::kInvalidIndex) {
                    context.log(4,"        Labeling neighbor %d with component %d and adding to queue.\n", v, comp_id);
                    final_component_label_[v] = comp_id; // Label neighbor
                    final_components_[comp_id].push_back(v); // Add to component list
                    cluster_queue_local.push_back(v); // Add to queue for exploration
                    // Check if neighbor is the designated root
                    if (v == context.root) {
                        root_component_index_ = comp_id;
                        context.log(4, "        Neighbor is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
                    }
                } else {
                     context.log(4,"        Neighbor %d already labeled or invalid index.\n", v);
                }
            }
        } else {
             context.log(4, "      Node %d index out of bounds for phase3_neighbors_.\n", u);
        }
    }
    context.log(4, "StrongPruner::label_component_recursive Exit for Component %d.\n", comp_id);
}

PCSTFast::IndexType StrongPruner::find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx) {
    context.log(3, "StrongPruner::find_best_component_root Entry: Component=%d\n", comp_idx);
    // Validate component index and check if component is empty
    if (static_cast<size_t>(comp_idx) >= final_components_.size() || final_components_[comp_idx].empty()) {
         context.log(1, "Warning: Invalid or empty component index %d in find_best_component_root.\n", comp_idx);
        return PCSTFast::kInvalidIndex;
    }
    const auto& comp_nodes = final_components_[comp_idx];
    PCSTFast::IndexType initial_root = comp_nodes[0]; // Pick the first node as an arbitrary initial root
    context.log(3, "  Using initial root %d for component %d.\n", initial_root, comp_idx);

    // Reset parent/payoff state specifically for nodes in this component
    context.log(4,"  Resetting strong pruning parent/payoff for component %d nodes.\n", comp_idx);
    size_t num_nodes = context.prizes.size();
    if(strong_pruning_parent_.size() != num_nodes) strong_pruning_parent_.resize(num_nodes);
    if(strong_pruning_payoff_.size() != num_nodes) strong_pruning_payoff_.resize(num_nodes);
    for (PCSTFast::IndexType node_idx : comp_nodes) {
        assert(static_cast<size_t>(node_idx) < num_nodes && "Node index out of bounds in component");
        strong_pruning_parent_[node_idx] = {PCSTFast::kInvalidIndex, 0.0};
        strong_pruning_payoff_[node_idx] = 0.0;
    }

    // 1. First DFS pass: Calculate payoffs rooted at the initial_root (don't mark deleted)
    context.log(4,"  Running initial strong_pruning_dfs from %d (MarkDeleted=false).\n", initial_root);
    strong_pruning_dfs(context, initial_root, false);

    // Check if initial root payoff calculation was successful
    if (static_cast<size_t>(initial_root) >= strong_pruning_payoff_.size()) {
          context.log(0,"Error: Initial root index %d out of bounds for payoff after DFS.\n", initial_root);
          assert(false && "Payoff index out of bounds after DFS");
          return PCSTFast::kInvalidIndex;
    }

    // Initialize best root found so far
    PCSTFast::IndexType best_root = initial_root;
    PCSTFast::ValueType best_value = strong_pruning_payoff_[initial_root];
    context.log(3, "  Initial Payoff at root %d = %.9g\n", best_root, best_value);

    // 2. Second pass: Propagate payoffs to find the node with the maximum payoff when acting as root
    context.log(4,"  Starting payoff propagation (rerooting) from initial root %d.\n", initial_root);
    propagate_payoffs_and_find_best(context, initial_root, best_root, best_value);

    context.log(3, "StrongPruner::find_best_component_root Exit. Best Root=%d, Best Payoff=%.9g\n", best_root, best_value);
    return best_root; // Return the index of the node found to have the best payoff
}

void StrongPruner::propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root_out, PCSTFast::ValueType& best_value_out) {
    context.log(4,"StrongPruner::propagate_payoffs Entry: InitialRoot=%d\n", initial_root);
    // Use a local stack for the traversal
    std::vector<PCSTFast::IndexType> strong_pruning_stack2_local;
    strong_pruning_stack2_local.reserve(context.prizes.size()); // Reserve capacity

    size_t num_nodes = context.prizes.size(); // Get total number of nodes

    // Add children of the initial root to start the propagation
    if (static_cast<size_t>(initial_root) < phase3_neighbors_.size()) {
         context.log(4, "  Adding children of initial root %d to propagation stack...\n", initial_root);
        for (const auto& edge : phase3_neighbors_[initial_root]) {
            PCSTFast::IndexType neighbor = edge.first;
            // Check if neighbor's parent is the initial root (meaning it's a child in the DFS tree)
            if (static_cast<size_t>(neighbor) < strong_pruning_parent_.size() && strong_pruning_parent_[neighbor].first == initial_root)
            {
                 context.log(4, "    Adding child %d to stack.\n", neighbor);
                strong_pruning_stack2_local.push_back(neighbor); // Add child to local stack
            }
        }
    }

    int nodes_processed = 0;
    // Process nodes using the stack (simulating recursion for rerooting)
    while (!strong_pruning_stack2_local.empty()) {
        PCSTFast::IndexType u = strong_pruning_stack2_local.back(); // Current node
        strong_pruning_stack2_local.pop_back();
        nodes_processed++;
        context.log(4,"    Processing node %d from propagation stack (item %d).\n", u, nodes_processed);

        // Basic validity checks
        if (static_cast<size_t>(u) >= num_nodes ||
            static_cast<size_t>(u) >= strong_pruning_parent_.size() ||
            static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
                 context.log(1,"Warning: Invalid index %d during payoff propagation.\n", u);
                 continue;
            }

        PCSTFast::IndexType p = strong_pruning_parent_[u].first; // Parent in the initial DFS tree
        PCSTFast::ValueType edge_cost = strong_pruning_parent_[u].second; // Cost of edge (u, p)
        context.log(4,"      Parent=%d, EdgeCost=%.9g\n", p, edge_cost);

        // Check parent validity
        if (p == PCSTFast::kInvalidIndex || static_cast<size_t>(p) >= strong_pruning_payoff_.size()) {
            context.log(1,"Warning: Invalid parent index %d for node %d during payoff propagation.\n", p, u);
            continue;
        }

        // --- Rerooting Calculation ---
        // payoff_u: Payoff of subtree rooted at u (calculated in first DFS pass relative to initial_root)
        PCSTFast::ValueType payoff_u = strong_pruning_payoff_[u];
        // payoff_p: Payoff of subtree rooted at p (calculated in first DFS pass relative to initial_root)
        PCSTFast::ValueType payoff_p = strong_pruning_payoff_[p];
        context.log(4,"      Payoffs: Subtree(Current)=%.9g, Parent(Rooted)=%.9g\n", payoff_u, payoff_p);

        // Contribution of u's subtree to p's payoff (only if positive net payoff)
        PCSTFast::ValueType u_contrib = (payoff_u > edge_cost) ? (payoff_u - edge_cost) : 0.0;
        context.log(4,"      Current contribution to parent = max(0, %.9g - %.9g) = %.9g\n", payoff_u, edge_cost, u_contrib);

        // Payoff of the tree excluding u's subtree (seen from p's perspective)
        PCSTFast::ValueType p_without_u = payoff_p - u_contrib;
        context.log(4,"      Parent payoff without current = %.9g - %.9g = %.9g\n", payoff_p, u_contrib, p_without_u);

        // Contribution of the rest of the tree (excluding u's subtree) to u's payoff if u were root
        PCSTFast::ValueType p_contrib = (p_without_u > edge_cost) ? (p_without_u - edge_cost) : 0.0;
        context.log(4,"      Parent contribution to current = max(0, %.9g - %.9g) = %.9g\n", p_without_u, edge_cost, p_contrib);

        // Total payoff if node u were the root of the component
        PCSTFast::ValueType u_total_payoff = payoff_u + p_contrib;
        context.log(4,"      Total payoff if node %d is root = %.9g + %.9g = %.9g\n", u, payoff_u, p_contrib, u_total_payoff);

        // Update best root found so far
        if (u_total_payoff > best_value_out) {
             context.log(3,"      New best root found: Node=%d, Payoff=%.9g (OldBest: Node=%d, Payoff=%.9g)\n", u, u_total_payoff, best_root_out, best_value_out);
            best_root_out = u;
            best_value_out = u_total_payoff;
        }

        // Update payoff for u to reflect its value if it were root (needed for children's calculations)
        strong_pruning_payoff_[u] = u_total_payoff;

        // Add children of current node u to the stack for further propagation
        if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
             context.log(4,"      Adding children of node %d to propagation stack...\n", u);
            for (const auto& edge : phase3_neighbors_[u]) {
                PCSTFast::IndexType v = edge.first; // Neighbor
                // Check if v is a child of u in the initial DFS tree
                if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u)
                {
                      context.log(4,"        Adding child %d to stack.\n", v);
                    strong_pruning_stack2_local.push_back(v); // Add child to local stack
                }
            }
        }
    }
    context.log(4, "StrongPruner::propagate_payoffs Exit. Processed %d nodes.\n", nodes_processed);
}

void StrongPruner::strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted) {
    context.log(3,"StrongPruner::strong_pruning_dfs Entry: Start=%d, MarkDeleted=%d\n", start_node, mark_deleted);
    // Use a local stack for the DFS traversal
    std::vector<std::pair<bool, PCSTFast::IndexType>> strong_pruning_stack_local;
    strong_pruning_stack_local.reserve(context.prizes.size()); // Reserve capacity

    size_t num_nodes = context.prizes.size(); // Get total number of nodes

    // Validate start node index
    if(static_cast<size_t>(start_node) >= num_nodes || static_cast<size_t>(start_node) >= strong_pruning_parent_.size()) {
         context.log(1,"Warning: Invalid start node %d for strong_pruning_dfs.\n", start_node);
         return;
    }

    // Initialize root of the DFS tree
    context.log(4,"  Initializing DFS root %d: Parent=Invalid, Pushing PreOrder.\n", start_node);
    strong_pruning_parent_[start_node] = {PCSTFast::kInvalidIndex, 0.0}; // Mark start node as root
    strong_pruning_stack_local.emplace_back(true, start_node); // Push {is_pre_order=true, node}

    int nodes_visited_pre = 0;
    int nodes_visited_post = 0;
    // Process stack until empty
    while(!strong_pruning_stack_local.empty()) {
        bool is_pre = strong_pruning_stack_local.back().first; // Is this a pre-order visit?
        PCSTFast::IndexType u = strong_pruning_stack_local.back().second; // Current node
        strong_pruning_stack_local.pop_back();

        // Check node validity
        if(static_cast<size_t>(u) >= num_nodes || static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
             context.log(1,"Warning: Invalid node index %d popped from DFS stack.\n", u);
             continue;
        }

        if(is_pre) { // Pre-order visit: Discover node, push children
            nodes_visited_pre++;
            context.log(4,"  DFS PreOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_pre);
            // Push post-order visit marker for this node
            context.log(4,"    Pushing PostOrder visit for node %d.\n", u);
            strong_pruning_stack_local.emplace_back(false, u);
            // Initialize payoff with own prize
            assert(static_cast<size_t>(u) < context.prizes.size() && "Node index out of bounds for prizes");
            strong_pruning_payoff_[u] = context.prizes[u];
            context.log(4,"    Initialized Payoff[%d] = Prize = %.9g\n", u, strong_pruning_payoff_[u]);

            // Explore neighbors using phase 3 adjacency list
            if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                 context.log(4,"    Exploring neighbors of %d...\n", u);
                for(const auto& edge : phase3_neighbors_[u]) {
                    PCSTFast::IndexType v = edge.first; // Neighbor node
                    context.log(4,"      Neighbor %d (Parent is %d).\n", v, strong_pruning_parent_[u].first);
                     // Don't go back up to the parent in the DFS tree
                     if(v == strong_pruning_parent_[u].first) {
                        context.log(4,"        Is parent, skipping.\n");
                        continue;
                     }
                     // Check neighbor validity
                     if(static_cast<size_t>(v) >= num_nodes) {
                          context.log(1,"Warning: Invalid neighbor index %d for node %d.\n", v, u);
                          continue;
                     }

                     // Mark neighbor's parent and push for pre-order visit
                      context.log(4,"      Setting Parent[%d] = %d (Cost=%.9g). Pushing PreOrder.\n", v, u, edge.second);
                     assert(static_cast<size_t>(v) < strong_pruning_parent_.size() && "Neighbor index out of bounds for parent array");
                     strong_pruning_parent_[v] = {u, edge.second};
                     strong_pruning_stack_local.emplace_back(true, v); // Push neighbor onto local stack
                }
            }
        } else { // Post-order visit: Process children's results, calculate final payoff for u
            nodes_visited_post++;
            context.log(4,"  DFS PostOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_post);
            // Aggregate results from children
            if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                 context.log(4,"    Aggregating child payoffs for node %d (CurrentPayoff=%.9g)\n", u, strong_pruning_payoff_[u]);
                for(const auto& edge : phase3_neighbors_[u]) {
                    PCSTFast::IndexType v = edge.first; // Neighbor
                    // Check if neighbor v is a child of u in the current DFS tree
                    if(static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u) {
                        context.log(4,"      Processing child %d...\n", v);
                        assert(static_cast<size_t>(v) < strong_pruning_payoff_.size() && "Child index out of bounds for payoff array");
                        PCSTFast::ValueType v_cost = strong_pruning_parent_[v].second; // Cost of edge (u,v)
                        // Calculate net payoff of the child's subtree (payoff - cost to connect)
                        PCSTFast::ValueType v_net_payoff = strong_pruning_payoff_[v] - v_cost;
                        context.log(4,"        ChildPayoff=%.9g, EdgeCost=%.9g -> NetChildPayoff=%.9g\n",
                                    strong_pruning_payoff_[v], v_cost, v_net_payoff);

                        // --- Pruning Condition ---
                        if(v_net_payoff <= 0.0) {
                             context.log(3,"        Child %d subtree payoff %.9g <= 0.\n", v, v_net_payoff);
                            // If payoff is non-positive, prune this subtree if requested
                            if(mark_deleted) {
                                 context.log(3, "        Strong Pruning: Pruning subtree at %d.\n", v);
                                // Mark node v and its descendants (excluding u) as deleted
                                mark_nodes_as_deleted_strong(context, v, u);
                            } else {
                                 context.log(4,"        (Not marking deleted as MarkDeleted=false)\n");
                            }
                            // Do not add non-positive payoff to the parent u
                        } else {
                            // If payoff is positive, add it to the parent's payoff
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

void StrongPruner::mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node) {
    context.log(4,"StrongPruner::mark_nodes_as_deleted_strong Entry: Start=%d, Parent=%d\n", start_node, parent_node);
    // Check if start node is already deleted or invalid
    if (static_cast<size_t>(start_node) >= node_deleted_.size() || node_deleted_[start_node]) {
        context.log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node);
        return;
    }

    // Use a local queue for the BFS traversal
    std::vector<PCSTFast::IndexType> cluster_queue_local;
    cluster_queue_local.reserve(context.prizes.size()); // Reserve capacity

    cluster_queue_local.push_back(start_node); // Start traversal
    node_deleted_[start_node] = 1; // Mark start node deleted (modifies member)
    context.log(3,"  Strong: Marking node %d and subtree (excluding parent %d) as deleted.\n", start_node, parent_node);

    size_t q_idx = 0;
    int count = 1;
    // Perform BFS
    while(q_idx < cluster_queue_local.size()) {
        PCSTFast::IndexType u = cluster_queue_local[q_idx++];
        context.log(4,"    Processing node %d from strong deletion queue (index %zu).\n", u, q_idx-1);
        // Explore neighbors in the phase 3 graph
        if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
            for(const auto& edge : phase3_neighbors_[u]) {
                PCSTFast::IndexType v = edge.first; // Neighbor
                // Don't traverse back to the parent node
                if(v == parent_node) continue;
                // If neighbor is valid and not deleted, mark it and add to queue
                if(static_cast<size_t>(v) < node_deleted_.size() && !node_deleted_[v]) {
                    node_deleted_[v] = 1; // Modify member
                    cluster_queue_local.push_back(v); // Use local queue
                    count++;
                    context.log(3,"    Strong: Deleted node %d (neighbor of %d).\n", v, u);
                }
            }
        }
    }
    context.log(4,"StrongPruner::mark_nodes_as_deleted_strong Exit. Deleted %d nodes.\n", count);
}

void StrongPruner::run_strong_pruning(const PruningContext& context) {
    context.log(3, "StrongPruner::run_strong_pruning Entry.\n");
    // Label connected components based on phase 3 graph
    label_final_components(context);
    context.log(2, "Strong Pruning: Found %zu components.\n", final_components_.size());

    // Process each component independently
    for (size_t comp_idx = 0; comp_idx < final_components_.size(); ++comp_idx) {
        if(final_components_[comp_idx].empty()) {
             context.log(3, "  Skipping empty component %zu.\n", comp_idx);
             continue;
        }
        context.log(2,"Strong Pruning: Processing component %zu (size %zu).\n", comp_idx, final_components_[comp_idx].size());
        PCSTFast::IndexType root_node;
        // Determine the root for this component's pruning DFS
        if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
            // If this component contains the designated root, use it
            root_node = context.root;
            context.log(3,"  Using designated root %d for component %zu.\n", root_node, comp_idx);
        } else {
            // Otherwise, find the node within the component that maximizes payoff
            root_node = find_best_component_root(context, static_cast<PCSTFast::IndexType>(comp_idx));
            context.log(3,"  Using best root %d for component %zu.\n", root_node, comp_idx);
        }

        // Run the pruning DFS from the chosen root if valid
        if(root_node != PCSTFast::kInvalidIndex && static_cast<size_t>(root_node) < context.prizes.size()) {
             context.log(3,"  Running final strong_pruning_dfs from %d (MarkDeleted=true).\n", root_node);
             // Reset parent/payoff for this specific DFS run (already done in find_best if not designated root)
             if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                   context.log(4,"  Resetting strong pruning parent/payoff for root component %zu nodes.\n", comp_idx);
                  for(PCSTFast::IndexType n : final_components_[comp_idx]) {
                       if(static_cast<size_t>(n) < strong_pruning_parent_.size()) strong_pruning_parent_[n] = {PCSTFast::kInvalidIndex, 0.0};
                       if(static_cast<size_t>(n) < strong_pruning_payoff_.size()) strong_pruning_payoff_[n] = 0.0;
                  }
             }
             strong_pruning_dfs(context, root_node, true); // Run DFS, marking nodes for deletion
        } else {
             context.log(1,"Warning: Skipping strong pruning for component %zu due to invalid root (%d).\n", comp_idx, root_node);
        }
    } // End loop through components

    // Filter phase 2 edges based on the final node_deleted_ flags
    phase3_result_local_.clear(); // Clear previous result
    phase3_result_local_.reserve(phase2_result_local_.size());
    context.log(2, "Strong pruning component processing complete. Filtering phase 2 edges based on deleted nodes...\n");
    for(PCSTFast::IndexType edge_idx : phase2_result_local_) { // Use IndexType
        // Check edge index validity
        if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;
        const auto& edge = context.edges[edge_idx];
        // Check node indices and deletion flags
        bool u_del = (static_cast<size_t>(edge.first) >= node_deleted_.size() || node_deleted_[edge.first]);
        bool v_del = (static_cast<size_t>(edge.second) >= node_deleted_.size() || node_deleted_[edge.second]);
        // Keep edge only if neither endpoint was deleted
        if (!u_del && !v_del) {
            phase3_result_local_.push_back(edge_idx);
        } else {
             context.log(4,"  Strong: Edge %d (%d,%d) removed due to deleted endpoint(s).\n", edge_idx, edge.first, edge.second);
        }
    }
    context.log(3, "StrongPruner::run_strong_pruning Exit. Final edge count %zu.\n", phase3_result_local_.size());
}

// Override setup to resize strong-specific members
void StrongPruner::setup(const PruningContext& context) {
    AdvancedPrunerBase::setup(context); // Call base setup first
    context.log(4, "StrongPruner::setup Entry (resizing strong-specific members).\n");
    // Resize members based on input size
    size_t num_nodes = context.prizes.size();
    final_component_label_.resize(num_nodes);
    final_components_.reserve(num_nodes); // Reserve, actual size determined by labeling
    strong_pruning_parent_.resize(num_nodes);
    strong_pruning_payoff_.resize(num_nodes);
     context.log(4, "StrongPruner::setup Exit.\n");
}

void StrongPruner::prune(const PruningContext& context,
           std::vector<PCSTFast::IndexType>& result_nodes,
           std::vector<PCSTFast::IndexType>& result_edges) {
    context.log(1, "Pruning: Strong. Setting up...\n");
    setup(context); // Includes base setup + strong-specific setup
    context.log(1, "Pruning: Running Strong pruning logic...\n");
    run_strong_pruning(context); // Run the specific strong pruning logic
    result_edges = phase3_result_local_; // Assign final Strong edges
    context.log(2, "Strong pruning complete. Building final node set...\n");
    build_pruned_node_set(context, result_nodes); // Build nodes based on final edges and deletion flags
    context.log(1, "Final Result (Strong Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
}


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
        case PCSTFast::PruningMethod::kUnknownPruning: // Explicitly handle unknown
        default:
            // Throw an exception for unsupported or unknown methods
            throw std::invalid_argument("Unsupported or unknown pruning method provided to factory.");
    }
}


} // namespace internal
} // namespace cluster_approx
