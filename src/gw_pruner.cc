#include "gw_pruner.h"
#include "pruning_context.h"

#include <vector>
#include <algorithm>
#include <cassert>

namespace cluster_approx {
    namespace internal {

        void GWPruner::mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index) {
            context.log(4, "GWPruner::mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);
            PCSTFast::IndexType current_cluster_index = start_node_index;
            int steps = 0;

            while (current_cluster_index != PCSTFast::kInvalidIndex && static_cast<size_t>(current_cluster_index) < context.clusters.size()) {
                steps++;
                if (static_cast<size_t>(current_cluster_index) >= cluster_necessary_local_.size()) {
                    context.log(2, "Warning: Cluster index %d out of bounds for necessary flag (%zu) during traversal.\n", current_cluster_index, cluster_necessary_local_.size());
                    break;
                }

                if (cluster_necessary_local_[current_cluster_index]) {
                    context.log(4, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
                    break;
                }

                context.log(4, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);
                cluster_necessary_local_[current_cluster_index] = true;
                current_cluster_index = context.clusters[current_cluster_index].merged_into;
                context.log(4, "  Moving up to parent cluster %d.\n", current_cluster_index);
            }
            context.log(4, "GWPruner::mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
        }

        void GWPruner::mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index) {
            context.log(4, "GWPruner::mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);

            if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
                context.log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
                return;
            }

            std::vector<PCSTFast::IndexType> cluster_queue_local;
            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node_index);
            node_deleted_[start_node_index] = 1;
            context.log(3, "  GW: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);

            size_t queue_index = 0;
            int nodes_deleted_count = 1;

            while (queue_index < cluster_queue_local.size()) {
                PCSTFast::IndexType current_node_index = cluster_queue_local[queue_index++];
                context.log(4,"    Processing node %d from GW deletion queue (index %zu).\n", current_node_index, queue_index-1);

                if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                     context.log(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
                    for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                        context.log(5,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);
                        if (neighbor_node_index == parent_node_index) {
                             context.log(5,"          Neighbor is parent, skipping.\n");
                            continue;
                        }

                        if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                            node_deleted_[neighbor_node_index] = 1;
                            cluster_queue_local.push_back(neighbor_node_index);
                            nodes_deleted_count++;
                            context.log(3, "    GW: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                        } else {
                             context.log(5,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
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
            phase3_result_local_.clear();
            phase3_result_local_.reserve(phase2_result_local_.size());
            cluster_necessary_local_.assign(context.clusters.size(), false);

            context.log(2, "Starting GW pruning reverse pass (processing %zu edges).\n", phase2_result_local_.size());
            int edges_kept = 0;
            int edges_discarded = 0;

            for (int ii = std::ssize(phase2_result_local_) - 1; ii >= 0; --ii) {
                PCSTFast::IndexType edge_idx = phase2_result_local_[ii];
                context.log(4, "GW reverse pass: Processing edge %d (index %d from end).\n", edge_idx, (int)phase2_result_local_.size()-1-ii);

                if(static_cast<size_t>(edge_idx) >= context.edges.size() ||
                   static_cast<size_t>(edge_idx) >= context.edge_info.size()) {
                    context.log(2,"Warning: Invalid edge index %d in GW prune loop.\n", edge_idx);
                    continue;
                }

                const auto& edge = context.edges[edge_idx];
                PCSTFast::IndexType uu = edge.first;
                PCSTFast::IndexType vv = edge.second;
                context.log(5, "  Edge %d connects nodes (%d, %d).\n", edge_idx, uu, vv);

                bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
                bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);
                 context.log(5, "  Node status: %d Deleted=%d, %d Deleted=%d.\n", uu, u_deleted, vv, v_deleted);

                if (u_deleted && v_deleted) {
                    context.log(3, "  GW: Both endpoints (%d, %d) of edge %d deleted. Skipping.\n", uu, vv, edge_idx);
                    edges_discarded++;
                    continue;
                }

                PCSTFast::IndexType inactive_merge_idx = context.edge_info[edge_idx].inactive_merge_event;
                context.log(5, "  Inactive merge event index for edge %d: %d.\n", edge_idx, inactive_merge_idx);

                if (inactive_merge_idx == PCSTFast::kInvalidIndex) {
                    context.log(3, "  GW: Edge %d (%d, %d) Active-Active. Keeping.\n", edge_idx, uu, vv);
                    phase3_result_local_.push_back(edge_idx);
                    edges_kept++;
                    context.log(4, "  Marking clusters necessary from nodes %d and %d (if not deleted).\n", uu, vv);
                    if (!u_deleted) mark_clusters_as_necessary_gw(context, uu);
                    if (!v_deleted) mark_clusters_as_necessary_gw(context, vv);
                }
                else {
                    if (static_cast<size_t>(inactive_merge_idx) >= context.inactive_merge_events.size()) {
                        context.log(0,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, edge_idx);
                        assert(false && "Invalid inactive merge index");
                        continue;
                    }
                    const auto& merge_event = context.inactive_merge_events[inactive_merge_idx];
                    PCSTFast::IndexType inactive_rep = merge_event.inactive_cluster_index;
                    PCSTFast::IndexType active_node = merge_event.active_cluster_node;
                    PCSTFast::IndexType inactive_node = merge_event.inactive_cluster_node;
                     context.log(4, "  Edge %d was Active-Inactive. ActiveNode=%d, InactiveNode=%d, InactiveRep=%d.\n",
                                edge_idx, active_node, inactive_node, inactive_rep);

                    bool is_inactive_necessary = (static_cast<size_t>(inactive_rep) < cluster_necessary_local_.size() && cluster_necessary_local_[inactive_rep]);
                    context.log(4, "  Checking necessity of inactive representative cluster %d: Necessary=%d\n",
                                inactive_rep, is_inactive_necessary);

                    if (is_inactive_necessary) {
                        context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d needed. Keeping.\n", edge_idx, uu, vv, inactive_rep);
                        phase3_result_local_.push_back(edge_idx);
                        edges_kept++;
                         context.log(4, "  Marking clusters necessary from nodes %d and %d.\n", active_node, inactive_node);
                        mark_clusters_as_necessary_gw(context, active_node);
                        mark_clusters_as_necessary_gw(context, inactive_node);
                    } else {
                        context.log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d not needed. Pruning from node %d.\n", edge_idx, uu, vv, inactive_rep, merge_event.inactive_cluster_node);
                        edges_discarded++;
                        mark_nodes_as_deleted_gw(context, inactive_node, active_node);
                    }
                }
            }

            context.log(2, "GW pruning reverse pass complete. Reversing phase 3 result.\n");
            std::reverse(phase3_result_local_.begin(), phase3_result_local_.end());
            context.log(3, "GWPruner::run_gw_pruning Exit. Final edge count: %zu (Kept: %d, Discarded: %d)\n", phase3_result_local_.size(), edges_kept, edges_discarded);
        }

        void GWPruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.log(3, "Pruning: GW. Setting up...\n");
            setup(context);
            cluster_necessary_local_.resize(context.clusters.size());
            context.log(3, "Pruning: Running GW pruning logic...\n");
            run_gw_pruning(context);
            result_edges = phase3_result_local_;
            context.log(2, "GW pruning complete. Building final node set...\n");
            build_pruned_node_set(context, result_nodes);
            context.log(3, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
