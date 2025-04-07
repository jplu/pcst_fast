#include "connect_final_pruner.h"
#include "pruning_context.h"
#include "logger.h"

#include <vector>
#include <set>
#include <queue>
#include <limits>
#include <algorithm>
#include <utility>
#include <cassert>
#include <tuple>

namespace cluster_approx {
    namespace internal {

        using internal::LogLevel;

        void ConnectFinalPruner::build_full_adjacency(const PruningContext& context,
            AdjacencyList& adj) const
        {
             context.logger.log(LogLevel::DEBUG, "ConnectFinalPruner::build_full_adjacency Entry.\n");
            size_t num_nodes = context.prizes.size();
            adj.assign(num_nodes, {});
            if (!context.edges.empty()) {
                size_t avg_degree_estimate = 2 * context.edges.size() / num_nodes + 1;
                for(auto& neighbors : adj) {
                    neighbors.reserve(avg_degree_estimate);
                }
            }


            for(size_t i=0; i < context.edges.size(); ++i) {
                 if (i >= context.costs.size()) {
                    context.logger.log(LogLevel::FATAL, "Error: Edge index %zu out of bounds for costs vector size %zu.\n", i, context.costs.size());
                    continue;
                 }
                const auto& edge = context.edges[i];
                if (static_cast<size_t>(edge.first) < num_nodes && static_cast<size_t>(edge.second) < num_nodes) {
                    PCSTFast::IndexType edge_idx = static_cast<PCSTFast::IndexType>(i);
                    adj[edge.first].emplace_back(edge.second, context.costs[i], edge_idx);
                    adj[edge.second].emplace_back(edge.first, context.costs[i], edge_idx);
                } else {
                    context.logger.log(LogLevel::WARNING, "Warning: Edge %zu connects out-of-bounds node (%d or %d).\n", i, edge.first, edge.second);
                }
            }
             context.logger.log(LogLevel::DEBUG, "ConnectFinalPruner::build_full_adjacency Exit.\n");
        }

        void ConnectFinalPruner::run_steiner_approximation(
            const PruningContext& context,
            const std::set<PCSTFast::IndexType>& target_nodes,
            const AdjacencyList& adj,
            std::set<PCSTFast::IndexType>& steiner_nodes_out,
            std::set<PCSTFast::IndexType>& steiner_edges_out
        ) {
             context.logger.log(LogLevel::INFO, "ConnectFinalPruner::run_steiner_approximation Entry. Targets=%zu\n", target_nodes.size());
             const size_t num_nodes = context.prizes.size();

            steiner_nodes_out.clear();
            steiner_edges_out.clear();

            if (target_nodes.empty()) {
                context.logger.log(LogLevel::INFO, "  No target nodes, returning empty Steiner tree.\n");
                return;
            }

            steiner_nodes_out = { *target_nodes.begin() };
            std::set<PCSTFast::IndexType> remaining_targets = target_nodes;
            remaining_targets.erase(remaining_targets.begin());
            context.logger.log(LogLevel::DEBUG, "  Starting Steiner approx with node %d. Remaining targets: %zu\n", *steiner_nodes_out.begin(), remaining_targets.size());

            using PQState = DijkstraState;
            std::priority_queue<PQState, std::vector<PQState>, std::greater<PQState>> pq;

            while (!remaining_targets.empty()) {
                 context.logger.log(LogLevel::DEBUG, "  Steiner Iteration: %zu targets remaining.\n", remaining_targets.size());
                 pq = {};
                std::vector<PCSTFast::ValueType> dist(num_nodes, std::numeric_limits<PCSTFast::ValueType>::infinity());
                std::vector<std::pair<PCSTFast::IndexType, PCSTFast::IndexType>> parent_info(
                    num_nodes, {PCSTFast::kInvalidIndex, PCSTFast::kInvalidIndex});

                 for (PCSTFast::IndexType tree_node : steiner_nodes_out) {
                    if (static_cast<size_t>(tree_node) >= num_nodes) {
                        context.logger.log(LogLevel::ERROR, "Error: Invalid node index %d from steiner_nodes_out during Dijkstra init.\n", tree_node);
                        continue;
                    }
                    dist[tree_node] = 0.0;
                    pq.push({0.0, tree_node, PCSTFast::kInvalidIndex});
                    context.logger.log(LogLevel::TRACE, "    Initializing PQ with node %d (dist 0).\n", tree_node);
                }

                PCSTFast::IndexType closest_target_found = PCSTFast::kInvalidIndex;
                PCSTFast::ValueType min_dist_to_target = std::numeric_limits<PCSTFast::ValueType>::infinity();

                while (!pq.empty()) {
                    PQState current = pq.top(); pq.pop();
                    PCSTFast::IndexType u = current.node;
                    PCSTFast::ValueType d = current.distance;
                    context.logger.log(LogLevel::TRACE, "    Dijkstra: Popped node %d (dist %.9g).\n", u, d);

                    if (static_cast<size_t>(u) >= num_nodes) {
                        context.logger.log(LogLevel::ERROR, "Error: Invalid node index %d popped from PQ.\n", u);
                        continue;
                    }

                    if (d > dist[u]) {
                         context.logger.log(LogLevel::TRACE, "      Skipping (already found shorter path %.9g).\n", dist[u]);
                        continue;
                    }

                    if (remaining_targets.count(u)) {
                        context.logger.log(LogLevel::DEBUG, "    Dijkstra: Reached remaining target node %d (dist %.9g).\n", u, d);
                        closest_target_found = u;
                        min_dist_to_target = d;
                        break;
                    }

                    if (closest_target_found != PCSTFast::kInvalidIndex && d >= min_dist_to_target) {
                             context.logger.log(LogLevel::TRACE, "      Pruning Dijkstra branch (dist %.9g >= current min %.9g).\n", d, min_dist_to_target);
                        continue;
                    }

                     if (static_cast<size_t>(u) < adj.size()) {
                         for (const auto& edge_tuple : adj[u]) {
                             PCSTFast::IndexType v = std::get<0>(edge_tuple);
                             PCSTFast::ValueType cost = std::get<1>(edge_tuple);
                             PCSTFast::IndexType edge_index = std::get<2>(edge_tuple);

                             if (static_cast<size_t>(v) >= num_nodes) {
                                  context.logger.log(LogLevel::ERROR, "Error: Invalid neighbor index %d for node %d via edge %d.\n", v, u, edge_index);
                                  continue;
                             }

                             PCSTFast::ValueType new_dist = d + cost;

                             if (new_dist < dist[v]) {
                                     context.logger.log(LogLevel::TRACE, "      Updating neighbor %d: new_dist=%.9g (via edge %d).\n", v, new_dist, edge_index);
                                     dist[v] = new_dist;
                                     parent_info[v] = {u, edge_index};
                                     pq.push({new_dist, v, edge_index});
                             }
                         }
                     }
                }

                if (closest_target_found != PCSTFast::kInvalidIndex) {
                    context.logger.log(LogLevel::INFO, "    Connecting closest target node %d (min dist %.9g) to Steiner tree.\n", closest_target_found, min_dist_to_target);
                    PCSTFast::IndexType curr = closest_target_found;
                    while (steiner_nodes_out.find(curr) == steiner_nodes_out.end()) {
                        if (static_cast<size_t>(curr) >= parent_info.size() || parent_info[curr].first == PCSTFast::kInvalidIndex) {
                             context.logger.log(LogLevel::ERROR,"Error: Path reconstruction failed for node %d. Parent info missing or invalid.\n", curr);
                             assert(false && "Path reconstruction failed");
                             closest_target_found = PCSTFast::kInvalidIndex;
                             break;
                        }

                        PCSTFast::IndexType parent_node = parent_info[curr].first;
                        PCSTFast::IndexType edge_to_add = parent_info[curr].second;
                        assert(edge_to_add != PCSTFast::kInvalidIndex && "Path reconstruction failed: edge index missing");

                        context.logger.log(LogLevel::DEBUG,"      Adding edge %d (parent %d -> curr %d) and node %d.\n", edge_to_add, parent_node, curr, curr);
                        steiner_edges_out.insert(edge_to_add);
                        steiner_nodes_out.insert(curr);

                        if (remaining_targets.count(curr)) {
                            context.logger.log(LogLevel::DEBUG,"        Node %d on path is also a target, removing from remaining.\n", curr);
                            remaining_targets.erase(curr);
                        }
                        curr = parent_node;
                    }

                    if (closest_target_found != PCSTFast::kInvalidIndex) {
                         context.logger.log(LogLevel::DEBUG,"      Path reconstruction complete. Connected to node %d in existing tree.\n", curr);
                         remaining_targets.erase(closest_target_found);
                    } else {
                         context.logger.log(LogLevel::ERROR, "Error: Steiner approximation failed during path reconstruction. Result may be incomplete.\n");
                         break;
                    }
                } else {
                     context.logger.log(LogLevel::WARNING,"Warning: Could not connect remaining targets (%zu). Steiner tree might be incomplete (graph disconnected?).\n", remaining_targets.size());
                    break;
                }
            }

             context.logger.log(LogLevel::INFO, "ConnectFinalPruner::run_steiner_approximation Exit. Nodes=%zu, Edges=%zu\n", steiner_nodes_out.size(), steiner_edges_out.size());
        }

        void ConnectFinalPruner::prune(
            const PruningContext& context,
            std::vector<PCSTFast::IndexType>& result_nodes,
            std::vector<PCSTFast::IndexType>& result_edges)
        {
             context.logger.log(LogLevel::INFO, "Pruning: ConnectFinalComponents (Steiner Approx).\n");
             result_nodes.clear();
             result_edges.clear();

            std::set<PCSTFast::IndexType> target_nodes_set;
            for(PCSTFast::IndexType i=0; i<static_cast<PCSTFast::IndexType>(context.prizes.size()); ++i) {
                bool is_good = (static_cast<size_t>(i) < context.node_good.size()) ? context.node_good[i] : false;
                 bool has_prize = context.prizes[i] > PCSTFast::kEpsilon;

                if (is_good || has_prize) {
                    target_nodes_set.insert(i);
                }
            }
            if (context.root != PCSTFast::kNoRoot) {
                 if (static_cast<size_t>(context.root) >= context.prizes.size()) {
                     context.logger.log(LogLevel::FATAL, "Error: Root node %d index out of bounds.\n", context.root);
                 } else {
                    target_nodes_set.insert(context.root);
                 }
            }

            if (target_nodes_set.empty()) {
                context.logger.log(LogLevel::WARNING,"No target (good/prize/root) nodes found. Returning empty result.\n");
                return;
            }
            context.logger.log(LogLevel::INFO, "Identified %zu target nodes to connect.\n", target_nodes_set.size());

            AdjacencyList full_adj;
            build_full_adjacency(context, full_adj);

            std::set<PCSTFast::IndexType> steiner_nodes_set;
            std::set<PCSTFast::IndexType> steiner_edges_set;
            run_steiner_approximation(context, target_nodes_set, full_adj, steiner_nodes_set, steiner_edges_set);

            result_edges.assign(steiner_edges_set.begin(), steiner_edges_set.end());
            result_nodes.assign(steiner_nodes_set.begin(), steiner_nodes_set.end());
            std::sort(result_nodes.begin(), result_nodes.end());
            std::sort(result_edges.begin(), result_edges.end());

             context.logger.log(LogLevel::INFO, "Final Result (ConnectFinal Steiner Approx): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
