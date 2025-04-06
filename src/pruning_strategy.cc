#include "pruning_strategy.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdarg>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <set>
#include <queue>
#include <map>
#include <limits>

namespace cluster_approx {
    namespace internal {
        void PruningContext::format_and_log(int level, const char* format, ...) const {
            if (logger && verbosity_level >= level) {
                constexpr int BUFFER_SIZE = 1024;
                char buffer[BUFFER_SIZE];
                va_list args;
                
                va_start(args, format);
                
                int needed = vsnprintf(buffer, BUFFER_SIZE, format, args);
                
                va_end(args);

                if (needed < 0) {
                    logger(level, "[Pruner Logging Error: vsnprintf failed]");
                } else if (static_cast<size_t>(needed) >= BUFFER_SIZE) {
                    buffer[BUFFER_SIZE - 1] = '\0';

                    if (BUFFER_SIZE > 4) { 
                        snprintf(buffer + BUFFER_SIZE - 5, 5, "..."); 
                    }

                    logger(level, std::string(buffer) + "(TRUNCATED)");
                } else {
                    logger(level, std::string(buffer));
                }
            }
        }

        void IPruner::build_node_set_base(const PruningContext& context, const std::vector<PCSTFast::IndexType>& edge_set, std::vector<PCSTFast::IndexType>& node_set) {
            context.format_and_log(4, "IPruner::build_node_set_base Entry (using %zu edges).\n", edge_set.size());
            node_set.clear();

            size_t num_nodes = context.prizes.size();

            node_set.reserve(num_nodes);

            std::vector<uint8_t> included_nodes_local(num_nodes, 0);

            for (PCSTFast::IndexType edge_idx : edge_set) {
                if (static_cast<size_t>(edge_idx) >= context.edges.size()) {
                    context.format_and_log(2, "Warning: Invalid edge index %d in build_node_set_base.\n", edge_idx);

                    continue;
                }

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

                if (is_good && !is_included) {
                    context.format_and_log(4, "  Adding isolated good node %d.\n", ii);
                    node_set.push_back(ii);
                }
            }

            std::sort(node_set.begin(), node_set.end());

            context.format_and_log(4, "IPruner::build_node_set_base Exit. Final node set size: %zu\n", node_set.size());
        }

        void NoPruner::prune(const PruningContext& context,
                    std::vector<PCSTFast::IndexType>& result_nodes,
                    std::vector<PCSTFast::IndexType>& result_edges) {
            context.format_and_log(3, "Pruning: None. Using Phase 1 result directly.\n");
            result_edges.assign(context.phase1_result.begin(), context.phase1_result.end());
            build_node_set_base(context, result_edges, result_nodes);
            context.format_and_log(3, "Final Result (No Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

        void SimplePruner::build_phase2(const PruningContext& context) {
            context.format_and_log(3, "SimplePruner::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
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
                        context.format_and_log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                    }
                } else {
                    context.format_and_log(2, "Warning: Invalid edge index %d in SimplePruner::build_phase2.\n", edge_idx);
                }
            }

            context.format_and_log(3, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
        }

        void SimplePruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.format_and_log(3, "Pruning: Simple. Running Phase 2 filtering.\n");
            build_phase2(context);

            result_edges = phase2_result_local_;

            build_node_set_base(context, result_edges, result_nodes);
            context.format_and_log(3, "Final Result (Simple Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

        void AdvancedPrunerBase::build_phase2(const PruningContext& context) {
            context.format_and_log(3, "AdvancedPrunerBase::build_phase2 Entry (Filtering %zu phase 1 edges).\n", context.phase1_result.size());
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
                        context.format_and_log(4, "  Phase 2 pruning: Removing edge %d (%d, %d) due to non-good endpoint(s).\n", edge_idx, edge.first, edge.second);
                    }
                } else {
                    context.format_and_log(2, "Warning: Invalid edge index %d in AdvancedPrunerBase::build_phase2.\n", edge_idx);
                }
            }

            context.format_and_log(3, "Pruning: Phase 2 (Connectivity). Edges remaining: %zu\n", phase2_result_local_.size());
        }

        void AdvancedPrunerBase::build_phase3_adjacency(const PruningContext& context) {
            context.format_and_log(3, "AdvancedPrunerBase::build_phase3_adjacency Entry (Using %zu phase 2 edges).\n", phase2_result_local_.size());

            size_t num_nodes = context.prizes.size();

            if (phase3_neighbors_.size() != num_nodes) {
                phase3_neighbors_.resize(num_nodes);
            }

            for (auto& neighbors : phase3_neighbors_) { neighbors.clear(); }

            int edges_added = 0;

            for (PCSTFast::IndexType edge_idx : phase2_result_local_) {
                if (static_cast<size_t>(edge_idx) < context.edges.size() && static_cast<size_t>(edge_idx) < context.costs.size()) {
                    const auto& edge = context.edges[edge_idx];
                    PCSTFast::ValueType cost = context.costs[edge_idx];

                    if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                        phase3_neighbors_[edge.first].emplace_back(edge.second, cost);
                        phase3_neighbors_[edge.second].emplace_back(edge.first, cost);

                        edges_added++;
                    } else {
                        context.format_and_log(0, "Error: Invalid node index in edge %d while building phase 3 adjacency list.\n", edge_idx);
                        assert(false && "Invalid node index in build_phase3_adjacency");
                    }
                } else {
                    context.format_and_log(2, "Warning: Invalid edge index %d found in phase2_result_ during adjacency build.\n", edge_idx);
                }
            }

            context.format_and_log(3, "AdvancedPrunerBase::build_phase3_adjacency Exit. Added %d edges (x2) to lists.\n", edges_added);
        }

        void AdvancedPrunerBase::build_pruned_node_set(const PruningContext& context, std::vector<PCSTFast::IndexType>& node_set) {
            context.format_and_log(3, "AdvancedPrunerBase::build_pruned_node_set Entry.\n");
            node_set.clear();

            size_t num_nodes = context.prizes.size();

            node_set.reserve(num_nodes);

            int nodes_included = 0;

            for (PCSTFast::IndexType ii = 0; ii < static_cast<PCSTFast::IndexType>(num_nodes); ++ii) {
                bool is_good = (static_cast<size_t>(ii) < context.node_good.size() && context.node_good[ii]);
                bool is_deleted = (static_cast<size_t>(ii) < node_deleted_.size() && node_deleted_[ii]);

                if (is_good && !is_deleted) {
                    node_set.push_back(ii);
                    nodes_included++;
                }
            }

            std::sort(node_set.begin(), node_set.end());
            context.format_and_log(3, "AdvancedPrunerBase::build_pruned_node_set Exit. Final node set size: %zu\n", static_cast<size_t>(nodes_included));
        }

        void AdvancedPrunerBase::setup(const PruningContext& context) {
            context.format_and_log(4, "AdvancedPrunerBase::setup Entry.\n");
            phase2_result_local_.reserve(context.phase1_result.size());

            size_t num_nodes = context.prizes.size();

            node_deleted_.assign(num_nodes, 0);
            phase3_neighbors_.resize(num_nodes);
            build_phase2(context);
            build_phase3_adjacency(context);
            context.format_and_log(4, "AdvancedPrunerBase::setup Exit.\n");
        }

        void GWPruner::mark_clusters_as_necessary_gw(const PruningContext& context, PCSTFast::IndexType start_node_index) {
            context.format_and_log(4, "GWPruner::mark_clusters_as_necessary_gw Entry: StartNode=%d\n", start_node_index);

            PCSTFast::IndexType current_cluster_index = start_node_index;
            int steps = 0;

            while (current_cluster_index != PCSTFast::kInvalidIndex && static_cast<size_t>(current_cluster_index) < context.clusters.size()) {
                steps++;

                if (static_cast<size_t>(current_cluster_index) >= cluster_necessary_local_.size()) {
                    context.format_and_log(2, "Warning: Cluster index %d out of bounds for necessary flag (%zu) during traversal.\n", current_cluster_index, cluster_necessary_local_.size());
                    break;
                }

                if (cluster_necessary_local_[current_cluster_index]) {
                    context.format_and_log(4, "  Cluster %d already marked necessary. Stopping traversal.\n", current_cluster_index);
                    break;
                }

                context.format_and_log(4, "  Marking cluster %d as necessary (Step %d).\n", current_cluster_index, steps);

                cluster_necessary_local_[current_cluster_index] = true;
                current_cluster_index = context.clusters[current_cluster_index].merged_into;

                context.format_and_log(4, "  Moving up to parent cluster %d.\n", current_cluster_index);
            }

            context.format_and_log(4, "GWPruner::mark_clusters_as_necessary_gw Exit (marked path from node %d).\n", start_node_index);
        }

        void GWPruner::mark_nodes_as_deleted_gw(const PruningContext& context, PCSTFast::IndexType start_node_index, PCSTFast::IndexType parent_node_index) {
            context.format_and_log(4, "GWPruner::mark_nodes_as_deleted_gw Entry: StartNode=%d, ParentNode=%d\n", start_node_index, parent_node_index);

            if (static_cast<size_t>(start_node_index) >= node_deleted_.size() || node_deleted_[start_node_index]) {
                context.format_and_log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node_index);

                return;
            }

            std::vector<PCSTFast::IndexType> cluster_queue_local;

            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node_index);

            node_deleted_[start_node_index] = 1;

            context.format_and_log(3, "  GW: Marking node %d and its subtree (excluding parent %d) as deleted.\n", start_node_index, parent_node_index);

            size_t queue_index = 0;
            int nodes_deleted_count = 1;

            while (queue_index < cluster_queue_local.size()) {
                PCSTFast::IndexType current_node_index = cluster_queue_local[queue_index++];

                context.format_and_log(4,"    Processing node %d from GW deletion queue (index %zu).\n", current_node_index, queue_index-1);

                if (static_cast<size_t>(current_node_index) < phase3_neighbors_.size()) {
                    context.format_and_log(4,"      Neighbors of %d: %zu\n", current_node_index, phase3_neighbors_[current_node_index].size());

                    for (const auto& [neighbor_node_index, cost] : phase3_neighbors_[current_node_index]) {
                        context.format_and_log(5,"        Checking neighbor %d (Parent is %d).\n", neighbor_node_index, parent_node_index);

                        if (neighbor_node_index == parent_node_index) {
                            context.format_and_log(5,"          Neighbor is parent, skipping.\n");

                            continue;
                        }

                        if (static_cast<size_t>(neighbor_node_index) < node_deleted_.size() && !node_deleted_[neighbor_node_index]) {
                            node_deleted_[neighbor_node_index] = 1;

                            cluster_queue_local.push_back(neighbor_node_index);

                            nodes_deleted_count++;

                            context.format_and_log(3, "    GW: Deleted node %d (neighbor of %d). Adding to queue.\n", neighbor_node_index, current_node_index);
                        } else {
                            context.format_and_log(5,"          Neighbor %d already deleted or invalid index.\n", neighbor_node_index);
                        }
                    }
                } else {
                    context.format_and_log(4, "      Node %d has no neighbors in phase 3 graph (or index out of bounds).\n", current_node_index);
                }
            }

            context.format_and_log(4, "GWPruner::mark_nodes_as_deleted_gw Exit. Deleted %d nodes starting from %d.\n", nodes_deleted_count, start_node_index);
        }

        void GWPruner::run_gw_pruning(const PruningContext& context) {
            context.format_and_log(3, "GWPruner::run_gw_pruning Entry.\n");
            phase3_result_local_.clear();
            phase3_result_local_.reserve(phase2_result_local_.size());
            cluster_necessary_local_.assign(context.clusters.size(), false);

            context.format_and_log(2, "Starting GW pruning reverse pass (processing %zu edges).\n", phase2_result_local_.size());

            int edges_kept = 0;
            int edges_discarded = 0;

            for (int ii = std::ssize(phase2_result_local_) - 1; ii >= 0; --ii) {
                PCSTFast::IndexType edge_idx = phase2_result_local_[ii];

                context.format_and_log(4, "GW reverse pass: Processing edge %d (index %d from end).\n", edge_idx, (int)phase2_result_local_.size()-1-ii);

                if(static_cast<size_t>(edge_idx) >= context.edges.size() ||

                static_cast<size_t>(edge_idx) >= context.edge_info.size()) {
                    context.format_and_log(2,"Warning: Invalid edge index %d in GW prune loop.\n", edge_idx);

                    continue;
                }

                const auto& edge = context.edges[edge_idx];
                PCSTFast::IndexType uu = edge.first;
                PCSTFast::IndexType vv = edge.second;

                context.format_and_log(5, "  Edge %d connects nodes (%d, %d).\n", edge_idx, uu, vv);

                bool u_deleted = (static_cast<size_t>(uu) >= node_deleted_.size() || node_deleted_[uu]);
                bool v_deleted = (static_cast<size_t>(vv) >= node_deleted_.size() || node_deleted_[vv]);

                context.format_and_log(5, "  Node status: %d Deleted=%d, %d Deleted=%d.\n", uu, u_deleted, vv, v_deleted);

                if (u_deleted && v_deleted) {
                    context.format_and_log(3, "  GW: Both endpoints (%d, %d) of edge %d deleted. Skipping.\n", uu, vv, edge_idx);

                    edges_discarded++;

                    continue;
                }

                PCSTFast::IndexType inactive_merge_idx = context.edge_info[edge_idx].inactive_merge_event;

                context.format_and_log(5, "  Inactive merge event index for edge %d: %d.\n", edge_idx, inactive_merge_idx);

                if (inactive_merge_idx == PCSTFast::kInvalidIndex) {
                    context.format_and_log(3, "  GW: Edge %d (%d, %d) Active-Active. Keeping.\n", edge_idx, uu, vv);
                    phase3_result_local_.push_back(edge_idx);

                    edges_kept++;

                    context.format_and_log(4, "  Marking clusters necessary from nodes %d and %d (if not deleted).\n", uu, vv);

                    if (!u_deleted) mark_clusters_as_necessary_gw(context, uu);
                    if (!v_deleted) mark_clusters_as_necessary_gw(context, vv);
                } else {
                    if (static_cast<size_t>(inactive_merge_idx) >= context.inactive_merge_events.size()) {
                        context.format_and_log(0,"Error: Invalid inactive merge event index %d for edge %d.\n", inactive_merge_idx, edge_idx);
                        assert(false && "Invalid inactive merge index");

                        continue;
                    }

                    const auto& merge_event = context.inactive_merge_events[inactive_merge_idx];
                    PCSTFast::IndexType inactive_rep = merge_event.inactive_cluster_index;
                    PCSTFast::IndexType active_node = merge_event.active_cluster_node;
                    PCSTFast::IndexType inactive_node = merge_event.inactive_cluster_node;

                    context.format_and_log(4, "  Edge %d was Active-Inactive. ActiveNode=%d, InactiveNode=%d, InactiveRep=%d.\n",
                                edge_idx, active_node, inactive_node, inactive_rep);

                    bool is_inactive_necessary = (static_cast<size_t>(inactive_rep) < cluster_necessary_local_.size() && cluster_necessary_local_[inactive_rep]);

                    context.format_and_log(4, "  Checking necessity of inactive representative cluster %d: Necessary=%d\n",
                                inactive_rep, is_inactive_necessary);

                    if (is_inactive_necessary) {
                        context.format_and_log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d needed. Keeping.\n", edge_idx, uu, vv, inactive_rep);
                        phase3_result_local_.push_back(edge_idx);

                        edges_kept++;

                        context.format_and_log(4, "  Marking clusters necessary from nodes %d and %d.\n", active_node, inactive_node);
                        mark_clusters_as_necessary_gw(context, active_node);
                        mark_clusters_as_necessary_gw(context, inactive_node);
                    } else {
                        context.format_and_log(3, "  GW: Edge %d (%d, %d) A-I, Inactive side C%d not needed. Pruning from node %d.\n", edge_idx, uu, vv, inactive_rep, merge_event.inactive_cluster_node);

                        edges_discarded++;

                        mark_nodes_as_deleted_gw(context, inactive_node, active_node);
                    }
                }
            }

            context.format_and_log(2, "GW pruning reverse pass complete. Reversing phase 3 result.\n");
            std::reverse(phase3_result_local_.begin(), phase3_result_local_.end());
            context.format_and_log(3, "GWPruner::run_gw_pruning Exit. Final edge count: %zu (Kept: %d, Discarded: %d)\n", phase3_result_local_.size(), edges_kept, edges_discarded);
        }

        void GWPruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.format_and_log(3, "Pruning: GW. Setting up...\n");
            setup(context);
            cluster_necessary_local_.resize(context.clusters.size());
            context.format_and_log(3, "Pruning: Running GW pruning logic...\n");
            run_gw_pruning(context);

            result_edges = phase3_result_local_;

            context.format_and_log(2, "GW pruning complete. Building final node set...\n");
            build_pruned_node_set(context, result_nodes);
            context.format_and_log(3, "Final Result (GW Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

        void StrongPruner::setup(const PruningContext& context) {
            AdvancedPrunerBase::setup(context);
            context.format_and_log(4, "StrongPruner::setup Entry (resizing strong-specific members).\n");

            size_t num_nodes = context.prizes.size();

            final_component_label_.resize(num_nodes);
            final_components_.reserve(num_nodes);
            strong_pruning_parent_.resize(num_nodes);
            strong_pruning_payoff_.resize(num_nodes);
            context.format_and_log(4, "StrongPruner::setup Exit.\n");
        }

        void StrongPruner::label_final_components(const PruningContext& context) {
            context.format_and_log(3, "StrongPruner::label_final_components Entry.\n");

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
                    PCSTFast::IndexType new_component_id = static_cast<PCSTFast::IndexType>(final_components_.size());
                    components_found++;

                    context.format_and_log(3, "  Found new component %d starting from node %d.\n", new_component_id, start_node);
                    final_components_.emplace_back();
                    label_component_bfs(context, start_node, new_component_id);
                    context.format_and_log(3, "  Component %d labeled. Size=%zu.\n", new_component_id, final_components_[new_component_id].size());
                }
            }

            context.format_and_log(3, "StrongPruner::label_final_components Exit. Found %d components.\n", components_found);
        }

        void StrongPruner::label_component_bfs(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id) {
            context.format_and_log(4, "StrongPruner::label_component_bfs Entry: Start=%d, ID=%d\n", start_node, comp_id);

            std::vector<PCSTFast::IndexType> cluster_queue_local;

            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node);
            assert(static_cast<size_t>(start_node) < final_component_label_.size() && "Start node index out of bounds for label");
            assert(static_cast<size_t>(comp_id) < final_components_.size() && "Component ID out of bounds for components vector");

            final_component_label_[start_node] = comp_id;

            final_components_[comp_id].push_back(start_node);
            context.format_and_log(4, "  Added start node %d to component %d.\n", start_node, comp_id);

            if (start_node == context.root) {
                root_component_index_ = comp_id;

                context.format_and_log(4, "  Start node is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
            }

            size_t q_idx = 0;

            while(q_idx < cluster_queue_local.size()) {
                PCSTFast::IndexType u = cluster_queue_local[q_idx++];

                context.format_and_log(4,"    Processing node %d from component labeling queue (index %zu).\n", u, q_idx-1);

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                    context.format_and_log(4,"      Neighbors of %d: %zu\n", u, phase3_neighbors_[u].size());

                    for(const auto& edge_pair : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge_pair.first;

                        if (static_cast<size_t>(v) < final_component_label_.size() && final_component_label_[v] == PCSTFast::kInvalidIndex) {
                            context.format_and_log(5,"        Labeling neighbor %d with component %d and adding to queue.\n", v, comp_id);

                            final_component_label_[v] = comp_id;

                            final_components_[comp_id].push_back(v);
                            cluster_queue_local.push_back(v);

                            if (v == context.root) {
                                root_component_index_ = comp_id;

                                context.format_and_log(5, "        Neighbor is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
                            }
                        } else {
                            context.format_and_log(5,"        Neighbor %d already labeled or invalid index.\n", v);
                        }
                    }
                } else {
                    context.format_and_log(4, "      Node %d index out of bounds for phase3_neighbors_.\n", u);
                }
            }

            context.format_and_log(4, "StrongPruner::label_component_bfs Exit for Component %d.\n", comp_id);
        }

        PCSTFast::IndexType StrongPruner::find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx) {
            context.format_and_log(3, "StrongPruner::find_best_component_root Entry: Component=%d\n", comp_idx);

            if (static_cast<size_t>(comp_idx) >= final_components_.size() || final_components_[comp_idx].empty()) {
                context.format_and_log(2, "Warning: Invalid or empty component index %d in find_best_component_root.\n", comp_idx);

                return PCSTFast::kInvalidIndex;
            }

            const auto& comp_nodes = final_components_[comp_idx];
            PCSTFast::IndexType initial_root = comp_nodes[0];

            context.format_and_log(3, "  Using initial root %d for component %d.\n", initial_root, comp_idx);
            context.format_and_log(4,"  Resetting strong pruning parent/payoff for component %d nodes.\n", comp_idx);

            size_t num_nodes = context.prizes.size();

            if(strong_pruning_parent_.size() != num_nodes) strong_pruning_parent_.resize(num_nodes);
            if(strong_pruning_payoff_.size() != num_nodes) strong_pruning_payoff_.resize(num_nodes);

            for (PCSTFast::IndexType node_idx : comp_nodes) {
                assert(static_cast<size_t>(node_idx) < num_nodes && "Node index out of bounds in component");

                strong_pruning_parent_[node_idx] = {PCSTFast::kInvalidIndex, 0.0};
                strong_pruning_payoff_[node_idx] = 0.0;
            }

            context.format_and_log(4,"  Running initial strong_pruning_dfs from %d (MarkDeleted=false).\n", initial_root);
            strong_pruning_dfs(context, initial_root, false);

            if (static_cast<size_t>(initial_root) >= strong_pruning_payoff_.size()) {
                context.format_and_log(0,"Error: Initial root index %d out of bounds for payoff after DFS.\n", initial_root);
                assert(false && "Payoff index out of bounds after DFS");

                return PCSTFast::kInvalidIndex;
            }

            PCSTFast::IndexType current_best_root = initial_root;
            PCSTFast::ValueType current_best_value = strong_pruning_payoff_[initial_root];

            context.format_and_log(3, "  Initial Payoff at root %d = %.9g\n", current_best_root, current_best_value);
            context.format_and_log(4,"  Starting payoff propagation (rerooting) from initial root %d.\n", initial_root);
            propagate_payoffs_and_find_best(context, initial_root, current_best_root, current_best_value);
            context.format_and_log(3, "StrongPruner::find_best_component_root Exit. Best Root=%d, Best Payoff=%.9g\n", current_best_root, current_best_value);

            return current_best_root;
        }

        void StrongPruner::propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root_out, PCSTFast::ValueType& best_value_out) {
            context.format_and_log(4,"StrongPruner::propagate_payoffs Entry: InitialRoot=%d\n", initial_root);

            std::vector<PCSTFast::IndexType> strong_pruning_stack2_local;

            strong_pruning_stack2_local.reserve(context.prizes.size());

            size_t num_nodes = context.prizes.size();

            if (static_cast<size_t>(initial_root) < phase3_neighbors_.size()) {
                context.format_and_log(4, "  Adding children of initial root %d to propagation stack...\n", initial_root);

                for (const auto& edge : phase3_neighbors_[initial_root]) {
                    PCSTFast::IndexType neighbor = edge.first;

                    if (static_cast<size_t>(neighbor) < strong_pruning_parent_.size() && strong_pruning_parent_[neighbor].first == initial_root)
                    {
                        context.format_and_log(5, "    Adding child %d to stack.\n", neighbor);
                        strong_pruning_stack2_local.push_back(neighbor);
                    }
                }
            }

            int nodes_processed = 0;

            while (!strong_pruning_stack2_local.empty()) {
                PCSTFast::IndexType u = strong_pruning_stack2_local.back();

                strong_pruning_stack2_local.pop_back();

                nodes_processed++;

                context.format_and_log(4,"    Processing node %d from propagation stack (item %d).\n", u, nodes_processed);

                if (static_cast<size_t>(u) >= num_nodes ||
                    static_cast<size_t>(u) >= strong_pruning_parent_.size() ||
                    static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
                        context.format_and_log(2,"Warning: Invalid index %d during payoff propagation.\n", u);

                        continue;
                    }

                PCSTFast::IndexType p = strong_pruning_parent_[u].first;
                PCSTFast::ValueType edge_cost = strong_pruning_parent_[u].second;

                context.format_and_log(5,"      Parent=%d, EdgeCost=%.9g\n", p, edge_cost);

                if (p == PCSTFast::kInvalidIndex || static_cast<size_t>(p) >= strong_pruning_payoff_.size()) {
                    context.format_and_log(2,"Warning: Invalid parent index %d for node %d during payoff propagation.\n", p, u);

                    continue;
                }

                PCSTFast::ValueType payoff_u = strong_pruning_payoff_[u];
                PCSTFast::ValueType payoff_p = strong_pruning_payoff_[p];

                context.format_and_log(5,"      Payoffs: Subtree(Current)=%.9g, Parent(Rooted)=%.9g\n", payoff_u, payoff_p);

                PCSTFast::ValueType u_contrib = (payoff_u > edge_cost) ? (payoff_u - edge_cost) : 0.0;

                context.format_and_log(5,"      Current contribution to parent = max(0, %.9g - %.9g) = %.9g\n", payoff_u, edge_cost, u_contrib);

                PCSTFast::ValueType p_without_u = payoff_p - u_contrib;

                context.format_and_log(5,"      Parent payoff without current = %.9g - %.9g = %.9g\n", payoff_p, u_contrib, p_without_u);

                PCSTFast::ValueType p_contrib = (p_without_u > edge_cost) ? (p_without_u - edge_cost) : 0.0;

                context.format_and_log(5,"      Parent contribution to current = max(0, %.9g - %.9g) = %.9g\n", p_without_u, edge_cost, p_contrib);

                PCSTFast::ValueType u_total_payoff = payoff_u + p_contrib;

                context.format_and_log(5,"      Total payoff if node %d is root = %.9g + %.9g = %.9g\n", u, payoff_u, p_contrib, u_total_payoff);

                if (u_total_payoff > best_value_out) {
                    context.format_and_log(3,"      New best root found: Node=%d, Payoff=%.9g (OldBest: Node=%d, Payoff=%.9g)\n", u, u_total_payoff, best_root_out, best_value_out);

                    best_root_out = u;
                    best_value_out = u_total_payoff;
                }

                strong_pruning_payoff_[u] = u_total_payoff;

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                    context.format_and_log(4,"      Adding children of node %d to propagation stack...\n", u);

                    for (const auto& edge : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge.first;

                        if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u)
                        {
                            context.format_and_log(5,"        Adding child %d to stack.\n", v);
                            strong_pruning_stack2_local.push_back(v);
                        }
                    }
                }
            }

            context.format_and_log(4, "StrongPruner::propagate_payoffs Exit. Processed %d nodes.\n", nodes_processed);
        }

        void StrongPruner::strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted) {
            context.format_and_log(3,"StrongPruner::strong_pruning_dfs Entry: Start=%d, MarkDeleted=%d\n", start_node, mark_deleted);

            std::vector<std::pair<bool, PCSTFast::IndexType>> strong_pruning_stack_local;

            strong_pruning_stack_local.reserve(context.prizes.size());

            size_t num_nodes = context.prizes.size();

            if (static_cast<size_t>(start_node) >= num_nodes ||
                static_cast<size_t>(start_node) >= strong_pruning_parent_.size()) {

                context.format_and_log(2,"Warning: Invalid start node %d for strong_pruning_dfs.\n", start_node);

                return;
            }

            context.format_and_log(4,"  Initializing DFS root %d: Parent=Invalid, Pushing PreOrder.\n", start_node);

            strong_pruning_parent_[start_node] = {PCSTFast::kInvalidIndex, 0.0};

            strong_pruning_stack_local.emplace_back(true, start_node);

            int nodes_visited_pre = 0;
            int nodes_visited_post = 0;

            while (!strong_pruning_stack_local.empty()) {
                bool is_pre = strong_pruning_stack_local.back().first;
                PCSTFast::IndexType u = strong_pruning_stack_local.back().second;

                strong_pruning_stack_local.pop_back();

                if (static_cast<size_t>(u) >= num_nodes || static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
                    context.format_and_log(2,"Warning: Invalid node index %d popped from DFS stack.\n", u);

                    continue;
                }

                if (is_pre) {
                    nodes_visited_pre++;

                    context.format_and_log(4,"  DFS PreOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_pre);
                    context.format_and_log(5,"    Pushing PostOrder visit for node %d.\n", u);
                    strong_pruning_stack_local.emplace_back(false, u);
                    assert(static_cast<size_t>(u) < context.prizes.size() && "Node index out of bounds for prizes");

                    strong_pruning_payoff_[u] = context.prizes[u];

                    context.format_and_log(5,"    Initialized Payoff[%d] = Prize = %.9g\n", u, strong_pruning_payoff_[u]);

                    if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                        context.format_and_log(4,"    Exploring neighbors of %d...\n", u);

                        for (const auto& edge : phase3_neighbors_[u]) {
                            PCSTFast::IndexType v = edge.first;

                            context.format_and_log(5,"      Neighbor %d (Parent is %d).\n", v, strong_pruning_parent_[u].first);

                            if (v == strong_pruning_parent_[u].first) {
                                context.format_and_log(5,"        Is parent, skipping.\n");

                                continue;
                            }

                            if (static_cast<size_t>(v) >= num_nodes) {
                                    context.format_and_log(2,"Warning: Invalid neighbor index %d for node %d.\n", v, u);

                                    continue;
                            }

                            context.format_and_log(5,"      Setting Parent[%d] = %d (Cost=%.9g). Pushing PreOrder.\n", v, u, edge.second);
                            assert(static_cast<size_t>(v) < strong_pruning_parent_.size() && "Neighbor index out of bounds for parent array");

                            strong_pruning_parent_[v] = {u, edge.second};

                            strong_pruning_stack_local.emplace_back(true, v);
                        }
                    }
                } else {
                    nodes_visited_post++;

                    context.format_and_log(4,"  DFS PostOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_post);

                    if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                        context.format_and_log(4,"    Aggregating child payoffs for node %d (CurrentPayoff=%.9g)\n", u, strong_pruning_payoff_[u]);

                        for (const auto& edge : phase3_neighbors_[u]) {
                            PCSTFast::IndexType v = edge.first;

                            if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u)
                            {
                                context.format_and_log(5,"      Processing child %d...\n", v);
                                assert(static_cast<size_t>(v) < strong_pruning_payoff_.size() && "Child index out of bounds for payoff array");

                                PCSTFast::ValueType v_cost = strong_pruning_parent_[v].second;
                                PCSTFast::ValueType v_net_payoff = strong_pruning_payoff_[v] - v_cost;

                                context.format_and_log(5,"        ChildPayoff=%.9g, EdgeCost=%.9g -> NetChildPayoff=%.9g\n",
                                            strong_pruning_payoff_[v], v_cost, v_net_payoff);

                                if (v_net_payoff <= 0.0) {
                                    context.format_and_log(3,"        Child %d subtree payoff %.9g <= 0.\n", v, v_net_payoff);

                                    if (mark_deleted) {
                                        context.format_and_log(3, "        Strong Pruning: Pruning subtree at %d.\n", v);
                                        mark_nodes_as_deleted_strong(context, v, u);
                                    } else {
                                        context.format_and_log(4,"        (Not marking deleted as MarkDeleted=false)\n");
                                    }
                                } else {
                                    context.format_and_log(5,"        Adding %.9g to Payoff[%d].\n", v_net_payoff, u);

                                    strong_pruning_payoff_[u] += v_net_payoff;
                                }
                            }
                        }

                        context.format_and_log(4,"    Final Payoff[%d]=%.9g after aggregating children.\n", u, strong_pruning_payoff_[u]);
                    }
                }
            }

            context.format_and_log(3, "StrongPruner::strong_pruning_dfs Exit. Visited %d nodes (pre), %d nodes (post).\n", nodes_visited_pre, nodes_visited_post);
        }

        void StrongPruner::mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node) {
            context.format_and_log(4,"StrongPruner::mark_nodes_as_deleted_strong Entry: Start=%d, Parent=%d\n", start_node, parent_node);

            if (static_cast<size_t>(start_node) >= node_deleted_.size() || node_deleted_[start_node]) {
                context.format_and_log(4, "  Node %d already deleted or invalid index. Returning.\n", start_node);

                return;
            }

            std::vector<PCSTFast::IndexType> cluster_queue_local;

            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node);

            node_deleted_[start_node] = 1;

            context.format_and_log(3,"  Strong: Marking node %d and subtree (excluding parent %d) as deleted.\n", start_node, parent_node);

            size_t q_idx = 0;
            int count = 1;

            while(q_idx < cluster_queue_local.size()) {
                PCSTFast::IndexType u = cluster_queue_local[q_idx++];

                context.format_and_log(4,"    Processing node %d from strong deletion queue (index %zu).\n", u, q_idx-1);

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                    context.format_and_log(4,"      Neighbors of %d: %zu\n", u, phase3_neighbors_[u].size());

                    for(const auto& edge : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge.first;

                        context.format_and_log(5,"        Checking neighbor %d (Parent is %d).\n", v, parent_node);

                        if(v == parent_node) { context.format_and_log(5,"          Is parent, skipping.\n"); continue; }

                        if(static_cast<size_t>(v) < node_deleted_.size() && !node_deleted_[v]) {
                            node_deleted_[v] = 1;

                            cluster_queue_local.push_back(v);

                            count++;

                            context.format_and_log(3,"    Strong: Deleted node %d (neighbor of %d).\n", v, u);
                        } else {
                            context.format_and_log(5,"          Neighbor %d already deleted or invalid index.\n", v);
                        }
                    }
                } else {
                    context.format_and_log(4, "      Node %d index out of bounds for phase3_neighbors_.\n", u);
                }
            }

            context.format_and_log(4,"StrongPruner::mark_nodes_as_deleted_strong Exit. Deleted %d nodes.\n", count);
        }

        void StrongPruner::run_strong_pruning(const PruningContext& context) {
            context.format_and_log(3, "StrongPruner::run_strong_pruning Entry.\n");
            label_final_components(context);
            context.format_and_log(2, "Strong Pruning: Found %zu components.\n", final_components_.size());

            for (size_t comp_idx = 0; comp_idx < final_components_.size(); ++comp_idx) {
                if(final_components_[comp_idx].empty()) { context.format_and_log(3, "  Skipping empty component %zu.\n", comp_idx); continue; }

                context.format_and_log(2,"Strong Pruning: Processing component %zu (size %zu).\n", comp_idx, final_components_[comp_idx].size());

                PCSTFast::IndexType root_node;

                if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                    root_node = context.root;

                    context.format_and_log(3,"  Using designated root %d for component %zu.\n", root_node, comp_idx);
                } else {
                    root_node = find_best_component_root(context, static_cast<PCSTFast::IndexType>(comp_idx));

                    context.format_and_log(3,"  Using best root %d for component %zu.\n", root_node, comp_idx);
                }

                if(root_node != PCSTFast::kInvalidIndex && static_cast<size_t>(root_node) < context.prizes.size()) {
                    context.format_and_log(3,"  Running final strong_pruning_dfs from %d (MarkDeleted=true).\n", root_node);

                    if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                            context.format_and_log(4,"  Resetting strong pruning parent/payoff for root component %zu nodes.\n", comp_idx);
                        for(PCSTFast::IndexType n : final_components_[comp_idx]) {
                                if(static_cast<size_t>(n) < strong_pruning_parent_.size()) strong_pruning_parent_[n] = {PCSTFast::kInvalidIndex, 0.0};
                                if(static_cast<size_t>(n) < strong_pruning_payoff_.size()) strong_pruning_payoff_[n] = 0.0;
                        }
                    }

                    strong_pruning_dfs(context, root_node, true);
                } else {
                    context.format_and_log(2,"Warning: Skipping strong pruning for component %zu due to invalid root (%d).\n", comp_idx, root_node);
                }
            }

            phase3_result_local_.clear();
            phase3_result_local_.reserve(phase2_result_local_.size());

            context.format_and_log(2, "Strong pruning component processing complete. Filtering phase 2 edges based on deleted nodes...\n");

            for(PCSTFast::IndexType edge_idx : phase2_result_local_) {
                if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;

                const auto& edge = context.edges[edge_idx];
                bool u_del = (static_cast<size_t>(edge.first) >= node_deleted_.size() || node_deleted_[edge.first]);
                bool v_del = (static_cast<size_t>(edge.second) >= node_deleted_.size() || node_deleted_[edge.second]);

                if (!u_del && !v_del) {
                    phase3_result_local_.push_back(edge_idx);
                } else {
                    context.format_and_log(4,"  Strong: Edge %d (%d,%d) removed due to deleted endpoint(s).\n", edge_idx, edge.first, edge.second);
                }
            }
            context.format_and_log(3, "StrongPruner::run_strong_pruning Exit. Final edge count %zu.\n", phase3_result_local_.size());
        }

        void StrongPruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.format_and_log(3, "Pruning: Strong. Setting up...\n");
            setup(context);
            context.format_and_log(3, "Pruning: Running Strong pruning logic...\n");
            run_strong_pruning(context);

            result_edges = phase3_result_local_;

            context.format_and_log(2, "Strong pruning complete. Building final node set...\n");
            build_pruned_node_set(context, result_nodes);
            context.format_and_log(3, "Final Result (Strong Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

        void ConnectFinalPruner::build_full_adjacency(const PruningContext& context,
            std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>>& adj) const
        {
            context.format_and_log(4, "ConnectFinalPruner::build_full_adjacency Entry.\n");

            size_t num_nodes = context.prizes.size();

            adj.assign(num_nodes, {});

            for(size_t i=0; i < context.edges.size(); ++i) {
                if (i >= context.costs.size()) {
                    context.format_and_log(0, "Error: Edge index %zu out of bounds for costs vector size %zu.\n", i, context.costs.size());

                    continue;
                }

                const auto& edge = context.edges[i];

                if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                    adj[edge.first].emplace_back(edge.second, context.costs[i]);
                    adj[edge.second].emplace_back(edge.first, context.costs[i]);
                } else {
                    context.format_and_log(2, "Warning: Edge %zu connects out-of-bounds node (%d or %d).\n", i, edge.first, edge.second);
                }
            }

            context.format_and_log(4, "ConnectFinalPruner::build_full_adjacency Exit.\n");
        }

        void ConnectFinalPruner::run_steiner_approximation(
            const PruningContext& context,
            const std::set<PCSTFast::IndexType>& target_nodes,
            const std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>>& adj,
            std::set<PCSTFast::IndexType>& steiner_nodes_out,
            std::set<PCSTFast::IndexType>& steiner_edges_out
        ) {
            context.format_and_log(3, "ConnectFinalPruner::run_steiner_approximation Entry. Targets=%zu\n", target_nodes.size());
            steiner_nodes_out.clear();
            steiner_edges_out.clear();

            if (target_nodes.empty()) {
                context.format_and_log(3, "  No target nodes, returning empty Steiner tree.\n");

                return;
            }

            steiner_nodes_out = { *target_nodes.begin() };
            std::set<PCSTFast::IndexType> remaining_targets = target_nodes;

            remaining_targets.erase(remaining_targets.begin());
            context.format_and_log(4, "  Starting Steiner approx with node %d. Remaining targets: %zu\n", *steiner_nodes_out.begin(), remaining_targets.size());

            using PQState = DijkstraState;
            std::priority_queue<PQState, std::vector<PQState>, std::greater<PQState>> pq;

            while (!remaining_targets.empty()) {
                context.format_and_log(4, "  Steiner Iteration: %zu targets remaining.\n", remaining_targets.size());

                pq = {};
                std::map<PCSTFast::IndexType, PCSTFast::ValueType> dist;
                std::map<PCSTFast::IndexType, std::pair<PCSTFast::IndexType, PCSTFast::IndexType>> parent_info;

                for (PCSTFast::IndexType tree_node : steiner_nodes_out) {
                    dist[tree_node] = 0.0;

                    pq.push({0.0, tree_node, PCSTFast::kInvalidIndex});

                    parent_info[tree_node] = {PCSTFast::kInvalidIndex, PCSTFast::kInvalidIndex};

                    context.format_and_log(5, "    Initializing PQ with node %d (dist 0).\n", tree_node);
                }

                PCSTFast::IndexType closest_target_found = PCSTFast::kInvalidIndex;
                PCSTFast::ValueType min_dist_to_target = std::numeric_limits<PCSTFast::ValueType>::infinity();

                while (!pq.empty()) {
                    PQState current = pq.top();

                    pq.pop();

                    PCSTFast::IndexType u = current.node;
                    PCSTFast::ValueType d = current.distance;

                    context.format_and_log(5, "    Dijkstra: Popped node %d (dist %.9g).\n", u, d);

                    if (dist.count(u) && d > dist[u]) {
                        context.format_and_log(5, "      Skipping (already found shorter path %.9g).\n", dist[u]);

                        continue;
                    }

                    if (remaining_targets.count(u)) {
                        context.format_and_log(4, "    Dijkstra: Reached remaining target node %d (dist %.9g).\n", u, d);

                        closest_target_found = u;
                        min_dist_to_target = d;

                        break;
                    }

                    if (closest_target_found != PCSTFast::kInvalidIndex && d >= min_dist_to_target) {
                            context.format_and_log(5, "      Pruning Dijkstra branch (dist %.9g >= current min %.9g).\n", d, min_dist_to_target);

                        continue;
                    }

                    if (static_cast<size_t>(u) < adj.size()) {
                        for (size_t i = 0; i < context.edges.size(); ++i) {
                            const auto& edge_pair = context.edges[i];
                            PCSTFast::IndexType v = PCSTFast::kInvalidIndex;

                            if (edge_pair.first == u) v = edge_pair.second;
                            else if (edge_pair.second == u) v = edge_pair.first;
                            else continue;

                            if (i >= context.costs.size()) continue;

                            PCSTFast::ValueType cost = context.costs[i];
                            PCSTFast::ValueType new_dist = d + cost;

                            if (!dist.count(v) || new_dist < dist[v]) {
                                    context.format_and_log(5, "      Updating neighbor %d: new_dist=%.9g (via edge %zu).\n", v, new_dist, i);

                                    dist[v] = new_dist;
                                    parent_info[v] = {u, static_cast<PCSTFast::IndexType>(i)};

                                    pq.push({new_dist, v, static_cast<PCSTFast::IndexType>(i)});
                            }
                        }
                    }
                }

                if (closest_target_found != PCSTFast::kInvalidIndex) {
                    context.format_and_log(3, "    Connecting closest target node %d (min dist %.9g) to Steiner tree.\n", closest_target_found, min_dist_to_target);

                    PCSTFast::IndexType curr = closest_target_found;

                    while (steiner_nodes_out.find(curr) == steiner_nodes_out.end()) {
                        assert(parent_info.count(curr) && "Path reconstruction failed: parent info missing");

                        PCSTFast::IndexType parent_node = parent_info[curr].first;
                        PCSTFast::IndexType edge_to_add = parent_info[curr].second;

                        assert(edge_to_add != PCSTFast::kInvalidIndex && "Path reconstruction failed: edge index missing");

                        context.format_and_log(4,"      Adding edge %d (parent %d -> curr %d) and node %d.\n", edge_to_add, parent_node, curr, curr);
                        steiner_edges_out.insert(edge_to_add);
                        steiner_nodes_out.insert(curr);

                        if (remaining_targets.count(curr)) {
                            context.format_and_log(4,"        Node %d on path is also a target, removing from remaining.\n", curr);
                            remaining_targets.erase(curr);
                        }

                        curr = parent_node;

                        assert(curr != PCSTFast::kInvalidIndex && "Path reconstruction failed: hit invalid parent");
                    }

                    context.format_and_log(4,"      Path reconstruction complete. Connected to node %d in existing tree.\n", curr);
                    remaining_targets.erase(closest_target_found);
                } else {
                    context.format_and_log(2,"Warning: Could not connect remaining targets (%zu). Steiner tree might be incomplete.\n", remaining_targets.size());

                    break;
                }
            }

            context.format_and_log(3, "ConnectFinalPruner::run_steiner_approximation Exit. Nodes=%zu, Edges=%zu\n", steiner_nodes_out.size(), steiner_edges_out.size());
        }

        void ConnectFinalPruner::prune(
            const PruningContext& context,
            std::vector<PCSTFast::IndexType>& result_nodes,
            std::vector<PCSTFast::IndexType>& result_edges)
        {
            context.format_and_log(3, "Pruning: ConnectFinalComponents (Steiner Approx).\n");
            result_nodes.clear();
            result_edges.clear();

            std::set<PCSTFast::IndexType> target_nodes_set;

            for(PCSTFast::IndexType i=0; i<static_cast<PCSTFast::IndexType>(context.prizes.size()); ++i) {
                bool is_good = (static_cast<size_t>(i) < context.node_good.size() && context.node_good[i]);
                bool has_prize = context.prizes[i] > PCSTFast::kEpsilon;

                if (is_good || has_prize) {
                    target_nodes_set.insert(i);
                }
            }

            if (context.root != PCSTFast::kNoRoot) {
                if (static_cast<size_t>(context.root) >= context.prizes.size()) {
                    context.format_and_log(0, "Error: Root node %d index out of bounds.\n", context.root);
                } else {
                    target_nodes_set.insert(context.root);
                }
            }

            if (target_nodes_set.empty()) {
                context.format_and_log(2,"No target (good/prize/root) nodes found. Returning empty result.\n");

                return;
            }

            context.format_and_log(2, "Identified %zu target nodes to connect.\n", target_nodes_set.size());

            std::vector<std::vector<std::pair<PCSTFast::IndexType, PCSTFast::ValueType>>> full_adj;

            build_full_adjacency(context, full_adj);

            std::set<PCSTFast::IndexType> steiner_nodes_set;
            std::set<PCSTFast::IndexType> steiner_edges_set;

            run_steiner_approximation(context, target_nodes_set, full_adj, steiner_nodes_set, steiner_edges_set);

            result_edges.assign(steiner_edges_set.begin(), steiner_edges_set.end());
            result_nodes.assign(steiner_nodes_set.begin(), steiner_nodes_set.end());
            std::sort(result_nodes.begin(), result_nodes.end());
            std::sort(result_edges.begin(), result_edges.end());

            context.format_and_log(3, "Final Result (ConnectFinal Steiner Approx): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

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
                case PCSTFast::PruningMethod::kConnectFinalComponents:
                    return std::make_unique<ConnectFinalPruner>();
                case PCSTFast::PruningMethod::kUnknownPruning:
                default:
                    throw std::invalid_argument("Unsupported or unknown pruning method provided to factory.");
            }
        }
    }
}
