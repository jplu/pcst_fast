// gw_pruner.cc
#include "gw_pruner.h"
#include "pruning_context.h"
#include "logger.h"

#include <vector>
#include <algorithm>
#include <cassert>

namespace cluster_approx {
    namespace internal {

        using internal::LogLevel;

        void GWPruner::mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index) {
            context.logger.log(LogLevel::DEBUG, "GWPruner::mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);
            PCSTFast::IndexType current_cluster_index = start_node_index;
            int steps = 0;

            while (current_cluster_index != PCSTFast::kInvalidIndex && static_cast<size_t>(current_cluster_index) < context.clusters.size()) {
                steps++;
                if (static_cast<size_t>(current_cluster_index) >= cluster_necessary_local_.size()) {
                    context.logger.log(LogLevel::WARNING, "Warning: Cluster index %d out of bounds for necessary flag (%zu) during traversal.\n", current_cluster_index, cluster_necessary_local_.size());
                    break;
                }

                if (cluster_necessary_local_[current_cluster_index]) {
                    context.logger.log(LogLevel::DEBUG, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
                    break;
                }

                context.logger.log(LogLevel::DEBUG, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);
                cluster_necessary_local_[current_cluster_index] = true;
                current_cluster_index = context.clusters[current_cluster_index].merged_into;
                context.logger.log(LogLevel::DEBUG, "  Moving up to parent cluster %d.\n", current_cluster_index);
            }
            context.logger.log(LogLevel::DEBUG, "GWPruner::mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
        }

        void GWPruner::mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index) {
            context.logger.log(LogLevel::DEBUG, "GWPruner::mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);

            if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
                context.logger.log(LogLevel::DEBUG, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);
                return;
            }

            std::vector<PCSTFast::IndexType> cluster_queue_local;
            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node_index);
            node_deleted_[start_node_index] = 1;
            context.logger.log(LogLevel::INFO, "  GW: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);

            size_t queue_index = 0;
            int nodes_deleted_count = 1;

            while (queue_index < cluster_queue_local.size()) {
                PCSTFast::IndexType current_node_index = cluster_queue_local[queue_index++];
                context.logger.log(LogLevel::DEBUG,"    Processing node %d from GW deletion queue (index %zu).\n", current_node_index, queue_index-1);

                if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                     context.logger.log(LogLevel::DEBUG,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());
                    for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                        context.logger.log(LogLevel::TRACE,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);
                        if (neighbor_node_index == parent_node_index) {
                             context.logger.log(LogLevel::TRACE,"          Neighbor is parent, skipping.\n");
                            continue;
                        }

                        if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                            node_deleted_[neighbor_node_index] = 1;
                            cluster_queue_local.push_back(neighbor_node_index);
                            nodes_deleted_count++;
                            context.logger.log(LogLevel::INFO, "    GW: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                        } else {
                             context.logger.log(LogLevel::TRACE,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
                        }
                    }
                } else {
                     context.logger.log(LogLevel::DEBUG, "      Node %d has no neighbors in phase 3 graph (or index out of bounds).\n", current_node_index);
                }
            }
            context.logger.log(LogLevel::DEBUG, "GWPruner::mark_nodes_as_deleted_gw Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
        }

        void GWPruner::run_gw_pruning(const PruningContext& context) {
            context.logger.log(LogLevel::INFO, "GWPruner::run_gw_pruning Entry.\n");
            phase3_result_local_.clear();
            phase3_result_local_.reserve(phase2_result_local_.size());
            cluster_necessary_local_.assign(context.clusters.size(), false);

            context.logger.log(LogLevel::WARNING, "Starting GW pruning reverse pass (processing %zu edges).\n", phase2_result_local_.size()); // Level 2 -> WARNING
            int edges_kept = 0;
            int edges_discarded = 0;

            for (int ii = std::ssize(phase2_result_local_) - 1; ii >= 0; --ii) {
                PCSTFast::IndexType edge_idx = phase2_result_local_[ii];
                context.logger.log(LogLevel::DEBUG, "GW reverse pass: Processing edge %d (index %d from end).\n", edge_idx, (int)phase2_result_local_.size()-1-ii);

                if(static_cast<size_t>(edge_idx) >= context.edges.size() ||
                   static_cast<size_t>(edge_idx) >= context.edge_info.size()) {
                    context.logger.log(LogLevel::WARNING,"Warning: Invalid edge index %d in GW prune loop.\n", edge_idx);
                    continue;
                }

                const auto& edge = context.edges[edge_idx];
                PCSTFast::IndexType uu = edge.first;
                PCSTFast::IndexType vv = edge.second;
                context.logger.log(LogLevel::TRACE, "  Edge %d connects nodes (%d, %d).\n", edge_idx, uu, vv);

                bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
                bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);
                 context.logger.log(LogLevel::TRACE, "  Node status: %d Deleted=%d, %d Deleted=%d.\n", uu, u_deleted, vv, v_deleted);

                if (u_deleted && v_deleted) {
                    context.logger.log(LogLevel::INFO, "  GW: Both endpoints (%d, %d) of edge %d deleted. Skipping.\n", uu, vv, edge_idx);
                    edges_discarded++;
                    continue;
                }

                PCSTFast::IndexType inactive_merge_idx = context.edge_info[edge_idx].inactive_merge_event;
                context.logger.log(LogLevel::TRACE, "  Inactive merge event index for edge %d: %d.\n", edge_idx, inactive_merge_idx);

                if (inactive_merge_idx == PCSTFast::kInvalidIndex) {
                    context.logger.log(LogLevel::INFO, "  GW: Edge %d (%d, %d) Active-Active. Keeping.\n", edge_idx, uu, vv);
                    phase3_result_local_.push_back(edge_idx);
                    edges_kept++;
                    context.logger.log(LogLevel::DEBUG, "  Marking clusters necessary from nodes %d and %d (if not deleted).\n", uu, vv);
                    if (!u_deleted) mark_clusters_as_necessary_gw(context, uu);
                    if (!v_deleted) mark_clusters_as_necessary_gw(context, vv);
                }
                else {
                    if (static_cast<size_t>(inactive_merge_idx) >= context.inactive_merge_events.size()) {
                        context.logger.log(LogLevel::FATAL,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, edge_idx);
                        assert(false && "Invalid inactive merge index");
                        continue;
                    }
                    const auto& merge_event = context.inactive_merge_events[inactive_merge_idx];
                    PCSTFast::IndexType inactive_rep = merge_event.inactive_cluster_index;
                    PCSTFast::IndexType active_node = merge_event.active_cluster_node;
                    PCSTFast::IndexType inactive_node = merge_event.inactive_cluster_node;
                     context.logger.log(LogLevel::DEBUG, "  Edge %d was Active-Inactive. ActiveNode=%d, InactiveNode=%d, InactiveRep=%d.\n",
                                edge_idx, active_node, inactive_node, inactive_rep);

                    bool is_inactive_necessary = (static_cast<size_t>(inactive_rep) < cluster_necessary_local_.size() && cluster_necessary_local_[inactive_rep]);
                    context.logger.log(LogLevel::DEBUG, "  Checking necessity of inactive representative cluster %d: Necessary=%d\n",
                                inactive_rep, is_inactive_necessary);

                    if (is_inactive_necessary) {
                        context.logger.log(LogLevel::INFO, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d needed. Keeping.\n", edge_idx, uu, vv, inactive_rep);
                        phase3_result_local_.push_back(edge_idx);
                        edges_kept++;
                         context.logger.log(LogLevel::DEBUG, "  Marking clusters necessary from nodes %d and %d.\n", active_node, inactive_node);
                        mark_clusters_as_necessary_gw(context, active_node);
                        mark_clusters_as_necessary_gw(context, inactive_node);
                    } else {
                        context.logger.log(LogLevel::INFO, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d not needed. Pruning from node %d.\n", edge_idx, uu, vv, inactive_rep, merge_event.inactive_cluster_node);
                        edges_discarded++;
                        mark_nodes_as_deleted_gw(context, inactive_node, active_node);
                    }
                }
            }

            context.logger.log(LogLevel::WARNING, "GW pruning reverse pass complete. Reversing phase 3 result.\n"); // Level 2 -> WARNING
            std::reverse(phase3_result_local_.begin(), phase3_result_local_.end());
            context.logger.log(LogLevel::INFO, "GWPruner::run_gw_pruning Exit. Final edge count: %zu (Kept: %d, Discarded: %d)\n", phase3_result_local_.size(), edges_kept, edges_discarded);
        }

        void GWPruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.logger.log(LogLevel::INFO, "Pruning: GW. Setting up...\n");
            setup(context);
            cluster_necessary_local_.resize(context.clusters.size());
            context.logger.log(LogLevel::INFO, "Pruning: Running GW pruning logic...\n");
            run_gw_pruning(context);
            result_edges = phase3_result_local_;
            context.logger.log(LogLevel::WARNING, "GW pruning complete. Building final node set...\n"); // Level 2 -> WARNING
            build_pruned_node_set(context, result_nodes);
            context.logger.log(LogLevel::INFO, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
