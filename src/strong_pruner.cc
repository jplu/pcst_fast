#include "strong_pruner.h"
#include "pruning_context.h"
#include "logger.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>

namespace cluster_approx {
    namespace internal {

        using internal::LogLevel;

        void StrongPruner::setup(const PruningContext& context) {
            AdvancedPrunerBase::setup(context);
            context.logger.log(LogLevel::DEBUG, "StrongPruner::setup Entry (resizing strong-specific members).\n");
            size_t num_nodes = context.prizes.size();
            final_component_label_.resize(num_nodes);
            final_components_.reserve(num_nodes);
            strong_pruning_parent_.resize(num_nodes);
            strong_pruning_payoff_.resize(num_nodes);
            context.logger.log(LogLevel::DEBUG, "StrongPruner::setup Exit.\n");
        }

        void StrongPruner::label_final_components(const PruningContext& context) {
            context.logger.log(LogLevel::INFO, "StrongPruner::label_final_components Entry.\n");
            size_t num_nodes = context.prizes.size();
            final_component_label_.assign(num_nodes, PCSTFast::kInvalidIndex);
            final_components_.clear();
            root_component_index_ = PCSTFast::kInvalidIndex;
            int components_found = 0;

            for (PCSTFast::IndexType start_node = 0; start_node < static_cast<PCSTFast::IndexType>(num_nodes); ++start_node) {
                bool has_neighbors = (static_cast<size_t>(start_node) < phase3_neighbors_.size() && !phase3_neighbors_[start_node].empty());
                bool is_isolated_good_node = (static_cast<size_t>(start_node) < context.node_good.size() && context.node_good[start_node])
                                           && (static_cast<size_t>(start_node) >= phase3_neighbors_.size() || phase3_neighbors_[start_node].empty());

                bool is_in_phase2_graph = has_neighbors || is_isolated_good_node;

                if (is_in_phase2_graph && static_cast<size_t>(start_node) < final_component_label_.size() && final_component_label_[start_node] == PCSTFast::kInvalidIndex) {
                    PCSTFast::IndexType new_component_id = static_cast<PCSTFast::IndexType>(final_components_.size());
                    components_found++;
                    context.logger.log(LogLevel::INFO, "  Found new component %d starting from node %d.\n", new_component_id, start_node);
                    final_components_.emplace_back();
                    label_component_bfs(context, start_node, new_component_id);
                    context.logger.log(LogLevel::INFO, "  Component %d labeled. Size=%zu.\n", new_component_id, final_components_[new_component_id].size());
                }
            }
            context.logger.log(LogLevel::INFO, "StrongPruner::label_final_components Exit. Found %d components.\n", components_found);
        }

        void StrongPruner::label_component_bfs(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType comp_id) {
             context.logger.log(LogLevel::DEBUG, "StrongPruner::label_component_bfs Entry: Start=%d, ID=%d\n", start_node, comp_id);
            std::vector<PCSTFast::IndexType> cluster_queue_local;
            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node);

            assert(static_cast<size_t>(start_node) < final_component_label_.size() && "Start node index out of bounds for label");
            assert(static_cast<size_t>(comp_id) < final_components_.size() && "Component ID out of bounds for components vector");

            final_component_label_[start_node] = comp_id;
            final_components_[comp_id].push_back(start_node);
             context.logger.log(LogLevel::DEBUG, "  Added start node %d to component %d.\n", start_node, comp_id);

            if (start_node == context.root) {
                root_component_index_ = comp_id;
                 context.logger.log(LogLevel::DEBUG, "  Start node is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
            }

            size_t q_idx = 0;
            while(q_idx < cluster_queue_local.size()) {
                PCSTFast::IndexType u = cluster_queue_local[q_idx++];
                 context.logger.log(LogLevel::DEBUG,"    Processing node %d from component labeling queue (index %zu).\n", u, q_idx-1);

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                     context.logger.log(LogLevel::DEBUG,"      Neighbors of %d: %zu\n", u, phase3_neighbors_[u].size());
                    for(const auto& edge_pair : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge_pair.first;
                        if (static_cast<size_t>(v) < final_component_label_.size() && final_component_label_[v] == PCSTFast::kInvalidIndex) {
                             context.logger.log(LogLevel::TRACE,"        Labeling neighbor %d with component %d and adding to queue.\n", v, comp_id);
                            final_component_label_[v] = comp_id;
                            final_components_[comp_id].push_back(v);
                            cluster_queue_local.push_back(v);

                            if (v == context.root) {
                                root_component_index_ = comp_id;
                                 context.logger.log(LogLevel::TRACE, "        Neighbor is designated root. Updated root_component_index_ = %d.\n", root_component_index_);
                            }
                        } else {
                              context.logger.log(LogLevel::TRACE,"        Neighbor %d already labeled or invalid index.\n", v);
                        }
                    }
                } else {
                     context.logger.log(LogLevel::DEBUG, "      Node %d index out of bounds for phase3_neighbors_.\n", u);
                }
            }
             context.logger.log(LogLevel::DEBUG, "StrongPruner::label_component_bfs Exit for Component %d.\n", comp_id);
        }

        PCSTFast::IndexType StrongPruner::find_best_component_root(const PruningContext& context, PCSTFast::IndexType comp_idx) {
            context.logger.log(LogLevel::INFO, "StrongPruner::find_best_component_root Entry: Component=%d\n", comp_idx);

            if (static_cast<size_t>(comp_idx) >= final_components_.size() || final_components_[comp_idx].empty()) {
                context.logger.log(LogLevel::WARNING, "Warning: Invalid or empty component index %d in find_best_component_root.\n", comp_idx);
                return PCSTFast::kInvalidIndex;
            }

            const auto& comp_nodes = final_components_[comp_idx];
            PCSTFast::IndexType initial_root = comp_nodes[0];
            context.logger.log(LogLevel::INFO, "  Using initial root %d for component %d.\n", initial_root, comp_idx);
            context.logger.log(LogLevel::DEBUG,"  Resetting strong pruning parent/payoff for component %d nodes.\n", comp_idx);

            size_t num_nodes = context.prizes.size();
            if(strong_pruning_parent_.size() != num_nodes) strong_pruning_parent_.resize(num_nodes);
            if(strong_pruning_payoff_.size() != num_nodes) strong_pruning_payoff_.resize(num_nodes);

            for (PCSTFast::IndexType node_idx : comp_nodes) {
                assert(static_cast<size_t>(node_idx) < num_nodes && "Node index out of bounds in component");
                strong_pruning_parent_[node_idx] = {PCSTFast::kInvalidIndex, 0.0};
                strong_pruning_payoff_[node_idx] = 0.0;
            }

             context.logger.log(LogLevel::DEBUG,"  Running initial strong_pruning_dfs from %d (MarkDeleted=false).\n", initial_root);
            strong_pruning_dfs(context, initial_root, false);

            if (static_cast<size_t>(initial_root) >= strong_pruning_payoff_.size()) {
                 context.logger.log(LogLevel::FATAL,"Error: Initial root index %d out of bounds for payoff after DFS.\n", initial_root);
                 assert(false && "Payoff index out of bounds after DFS");
                return PCSTFast::kInvalidIndex;
            }

            PCSTFast::IndexType current_best_root = initial_root;
            PCSTFast::ValueType current_best_value = strong_pruning_payoff_[initial_root];
             context.logger.log(LogLevel::INFO, "  Initial Payoff at root %d = %.9g\n", current_best_root, current_best_value);

             context.logger.log(LogLevel::DEBUG,"  Starting payoff propagation (rerooting) from initial root %d.\n", initial_root);
            propagate_payoffs_and_find_best(context, initial_root, current_best_root, current_best_value);
            context.logger.log(LogLevel::INFO, "StrongPruner::find_best_component_root Exit. Best Root=%d, Best Payoff=%.9g\n", current_best_root, current_best_value);
            return current_best_root;
        }

        void StrongPruner::propagate_payoffs_and_find_best(const PruningContext& context, PCSTFast::IndexType initial_root, PCSTFast::IndexType& best_root_out, PCSTFast::ValueType& best_value_out) {
            context.logger.log(LogLevel::DEBUG,"StrongPruner::propagate_payoffs Entry: InitialRoot=%d\n", initial_root);
            std::vector<PCSTFast::IndexType> strong_pruning_stack2_local;
            strong_pruning_stack2_local.reserve(context.prizes.size());

            size_t num_nodes = context.prizes.size();
            if (static_cast<size_t>(initial_root) < phase3_neighbors_.size()) {
                context.logger.log(LogLevel::DEBUG, "  Adding children of initial root %d to propagation stack...\n", initial_root);
                for (const auto& edge : phase3_neighbors_[initial_root]) {
                    PCSTFast::IndexType neighbor = edge.first;
                    if (static_cast<size_t>(neighbor) < strong_pruning_parent_.size() && strong_pruning_parent_[neighbor].first == initial_root)
                    {
                         context.logger.log(LogLevel::TRACE, "    Adding child %d to stack.\n", neighbor);
                        strong_pruning_stack2_local.push_back(neighbor);
                    }
                }
            }

            int nodes_processed = 0;
            while (!strong_pruning_stack2_local.empty()) {
                PCSTFast::IndexType u = strong_pruning_stack2_local.back();
                strong_pruning_stack2_local.pop_back();
                nodes_processed++;
                 context.logger.log(LogLevel::DEBUG,"    Processing node %d from propagation stack (item %d).\n", u, nodes_processed);

                 if (static_cast<size_t>(u) >= num_nodes ||
                     static_cast<size_t>(u) >= strong_pruning_parent_.size() ||
                     static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
                         context.logger.log(LogLevel::WARNING,"Warning: Invalid index %d during payoff propagation.\n", u);
                         continue;
                     }

                PCSTFast::IndexType p = strong_pruning_parent_[u].first;
                PCSTFast::ValueType edge_cost = strong_pruning_parent_[u].second;
                 context.logger.log(LogLevel::TRACE,"      Parent=%d, EdgeCost=%.9g\n", p, edge_cost);

                 if (p == PCSTFast::kInvalidIndex || static_cast<size_t>(p) >= strong_pruning_payoff_.size()) {
                     context.logger.log(LogLevel::WARNING,"Warning: Invalid parent index %d for node %d during payoff propagation.\n", p, u);
                     continue;
                 }

                PCSTFast::ValueType payoff_u = strong_pruning_payoff_[u];
                PCSTFast::ValueType payoff_p = strong_pruning_payoff_[p];
                 context.logger.log(LogLevel::TRACE,"      Payoffs: Subtree(Current)=%.9g, Parent(Rooted)=%.9g\n", payoff_u, payoff_p);

                PCSTFast::ValueType u_contrib = (payoff_u > edge_cost) ? (payoff_u - edge_cost) : 0.0;
                 context.logger.log(LogLevel::TRACE,"      Current contribution to parent = max(0, %.9g - %.9g) = %.9g\n", payoff_u, edge_cost, u_contrib);

                PCSTFast::ValueType p_without_u = payoff_p - u_contrib;
                 context.logger.log(LogLevel::TRACE,"      Parent payoff without current = %.9g - %.9g = %.9g\n", payoff_p, u_contrib, p_without_u);

                PCSTFast::ValueType p_contrib = (p_without_u > edge_cost) ? (p_without_u - edge_cost) : 0.0;
                 context.logger.log(LogLevel::TRACE,"      Parent contribution to current = max(0, %.9g - %.9g) = %.9g\n", p_without_u, edge_cost, p_contrib);

                PCSTFast::ValueType u_total_payoff = payoff_u + p_contrib;
                 context.logger.log(LogLevel::TRACE,"      Total payoff if node %d is root = %.9g + %.9g = %.9g\n", u, payoff_u, p_contrib, u_total_payoff);

                if (u_total_payoff > best_value_out) {
                     context.logger.log(LogLevel::INFO,"      New best root found: Node=%d, Payoff=%.9g (OldBest: Node=%d, Payoff=%.9g)\n", u, u_total_payoff, best_root_out, best_value_out);
                    best_root_out = u;
                    best_value_out = u_total_payoff;
                }

                strong_pruning_payoff_[u] = u_total_payoff;

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                     context.logger.log(LogLevel::DEBUG,"      Adding children of node %d to propagation stack...\n", u);
                    for (const auto& edge : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge.first;
                        if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u)
                        {
                             context.logger.log(LogLevel::TRACE,"        Adding child %d to stack.\n", v);
                            strong_pruning_stack2_local.push_back(v);
                        }
                    }
                }
            }
             context.logger.log(LogLevel::DEBUG, "StrongPruner::propagate_payoffs Exit. Processed %d nodes.\n", nodes_processed);
        }

        void StrongPruner::strong_pruning_dfs(const PruningContext& context, PCSTFast::IndexType start_node, bool mark_deleted) {
             context.logger.log(LogLevel::INFO,"StrongPruner::strong_pruning_dfs Entry: Start=%d, MarkDeleted=%d\n", start_node, mark_deleted);
            std::vector<std::pair<bool, PCSTFast::IndexType>> strong_pruning_stack_local;
            strong_pruning_stack_local.reserve(context.prizes.size());

            size_t num_nodes = context.prizes.size();
             if (static_cast<size_t>(start_node) >= num_nodes ||
                 static_cast<size_t>(start_node) >= strong_pruning_parent_.size()) {
                 context.logger.log(LogLevel::WARNING,"Warning: Invalid start node %d for strong_pruning_dfs.\n", start_node);
                 return;
             }

             context.logger.log(LogLevel::DEBUG,"  Initializing DFS root %d: Parent=Invalid, Pushing PreOrder.\n", start_node);
            strong_pruning_parent_[start_node] = {PCSTFast::kInvalidIndex, 0.0};
            strong_pruning_stack_local.emplace_back(true, start_node);

            int nodes_visited_pre = 0;
            int nodes_visited_post = 0;

            while (!strong_pruning_stack_local.empty()) {
                bool is_pre = strong_pruning_stack_local.back().first;
                PCSTFast::IndexType u = strong_pruning_stack_local.back().second;
                strong_pruning_stack_local.pop_back();

                if (static_cast<size_t>(u) >= num_nodes || static_cast<size_t>(u) >= strong_pruning_payoff_.size()) {
                    context.logger.log(LogLevel::WARNING,"Warning: Invalid node index %d popped from DFS stack.\n", u);
                    continue;
                }

                if (is_pre) {
                    nodes_visited_pre++;
                     context.logger.log(LogLevel::DEBUG,"  DFS PreOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_pre);
                     context.logger.log(LogLevel::TRACE,"    Pushing PostOrder visit for node %d.\n", u);
                    strong_pruning_stack_local.emplace_back(false, u);
                     assert(static_cast<size_t>(u) < context.prizes.size() && "Node index out of bounds for prizes");
                    strong_pruning_payoff_[u] = context.prizes[u];
                     context.logger.log(LogLevel::TRACE,"    Initialized Payoff[%d] = Prize = %.9g\n", u, strong_pruning_payoff_[u]);

                    if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                         context.logger.log(LogLevel::DEBUG,"    Exploring neighbors of %d...\n", u);
                        for (const auto& edge : phase3_neighbors_[u]) {
                            PCSTFast::IndexType v = edge.first;
                             context.logger.log(LogLevel::TRACE,"      Neighbor %d (Parent is %d).\n", v, strong_pruning_parent_[u].first);

                            if (v == strong_pruning_parent_[u].first) {
                                 context.logger.log(LogLevel::TRACE,"        Is parent, skipping.\n");
                                continue;
                            }
                             if (static_cast<size_t>(v) >= num_nodes) {
                                     context.logger.log(LogLevel::WARNING,"Warning: Invalid neighbor index %d for node %d.\n", v, u);
                                     continue;
                             }

                             context.logger.log(LogLevel::TRACE,"      Setting Parent[%d] = %d (Cost=%.9g). Pushing PreOrder.\n", v, u, edge.second);
                             assert(static_cast<size_t>(v) < strong_pruning_parent_.size() && "Neighbor index out of bounds for parent array");
                            strong_pruning_parent_[v] = {u, edge.second};
                            strong_pruning_stack_local.emplace_back(true, v);
                        }
                    }
                }
                else {
                    nodes_visited_post++;
                    context.logger.log(LogLevel::DEBUG,"  DFS PostOrder Visit: Node=%d (Visit %d)\n", u, nodes_visited_post);

                    if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                        context.logger.log(LogLevel::DEBUG,"    Aggregating child payoffs for node %d (CurrentPayoff=%.9g)\n", u, strong_pruning_payoff_[u]);
                        for (const auto& edge : phase3_neighbors_[u]) {
                            PCSTFast::IndexType v = edge.first;
                            if (static_cast<size_t>(v) < strong_pruning_parent_.size() && strong_pruning_parent_[v].first == u)
                            {
                                 context.logger.log(LogLevel::TRACE,"      Processing child %d...\n", v);
                                 assert(static_cast<size_t>(v) < strong_pruning_payoff_.size() && "Child index out of bounds for payoff array");
                                PCSTFast::ValueType v_cost = strong_pruning_parent_[v].second;
                                PCSTFast::ValueType v_net_payoff = strong_pruning_payoff_[v] - v_cost;
                                 context.logger.log(LogLevel::TRACE,"        ChildPayoff=%.9g, EdgeCost=%.9g -> NetChildPayoff=%.9g\n",
                                             strong_pruning_payoff_[v], v_cost, v_net_payoff);

                                if (v_net_payoff <= 0.0) {
                                     context.logger.log(LogLevel::INFO,"        Child %d subtree payoff %.9g <= 0.\n", v, v_net_payoff);
                                    if (mark_deleted) {
                                         context.logger.log(LogLevel::INFO, "        Strong Pruning: Pruning subtree at %d.\n", v);
                                        mark_nodes_as_deleted_strong(context, v, u);
                                    } else {
                                         context.logger.log(LogLevel::DEBUG,"        (Not marking deleted as MarkDeleted=false)\n");
                                    }
                                } else {
                                     context.logger.log(LogLevel::TRACE,"        Adding %.9g to Payoff[%d].\n", v_net_payoff, u);
                                    strong_pruning_payoff_[u] += v_net_payoff;
                                }
                            }
                        }
                         context.logger.log(LogLevel::DEBUG,"    Final Payoff[%d]=%.9g after aggregating children.\n", u, strong_pruning_payoff_[u]);
                    }
                }
            }
             context.logger.log(LogLevel::INFO, "StrongPruner::strong_pruning_dfs Exit. Visited %d nodes (pre), %d nodes (post).\n", nodes_visited_pre, nodes_visited_post);
        }

        void StrongPruner::mark_nodes_as_deleted_strong(const PruningContext& context, PCSTFast::IndexType start_node, PCSTFast::IndexType parent_node) {
            context.logger.log(LogLevel::DEBUG,"StrongPruner::mark_nodes_as_deleted_strong Entry: Start=%d, Parent=%d\n", start_node, parent_node);
            if (static_cast<size_t>(start_node) >= node_deleted_.size() || node_deleted_[start_node]) {
                context.logger.log(LogLevel::DEBUG, "  Node %d already deleted or invalid index. Returning.\n", start_node);
                return;
            }

            std::vector<PCSTFast::IndexType> cluster_queue_local;
            cluster_queue_local.reserve(context.prizes.size());
            cluster_queue_local.push_back(start_node);
            node_deleted_[start_node] = 1;
             context.logger.log(LogLevel::INFO,"  Strong: Marking node %d and subtree (excluding parent %d) as deleted.\n", start_node, parent_node);

            size_t q_idx = 0;
            int count = 1;
            while(q_idx < cluster_queue_local.size()) {
                PCSTFast::IndexType u = cluster_queue_local[q_idx++];
                 context.logger.log(LogLevel::DEBUG,"    Processing node %d from strong deletion queue (index %zu).\n", u, q_idx-1);

                if (static_cast<size_t>(u) < phase3_neighbors_.size()) {
                     context.logger.log(LogLevel::DEBUG,"      Neighbors of %d: %zu\n", u, phase3_neighbors_[u].size());
                    for(const auto& edge : phase3_neighbors_[u]) {
                        PCSTFast::IndexType v = edge.first;
                        context.logger.log(LogLevel::TRACE,"        Checking neighbor %d (Parent is %d).\n", v, parent_node);
                        if(v == parent_node) { context.logger.log(LogLevel::TRACE,"          Is parent, skipping.\n"); continue; }
                        if(static_cast<size_t>(v) < node_deleted_.size() && !node_deleted_[v]) {
                            node_deleted_[v] = 1;
                            cluster_queue_local.push_back(v);
                            count++;
                            context.logger.log(LogLevel::INFO,"    Strong: Deleted node %d (neighbor of %d).\n", v, u);
                        } else {
                            context.logger.log(LogLevel::TRACE,"          Neighbor %d already deleted or invalid index.\n", v);
                        }
                    }
                } else {
                     context.logger.log(LogLevel::DEBUG, "      Node %d index out of bounds for phase3_neighbors_.\n", u);
                }
            }
            context.logger.log(LogLevel::DEBUG,"StrongPruner::mark_nodes_as_deleted_strong Exit. Deleted %d nodes.\n", count);
        }

        void StrongPruner::run_strong_pruning(const PruningContext& context) {
            context.logger.log(LogLevel::INFO, "StrongPruner::run_strong_pruning Entry.\n");
            label_final_components(context);
            context.logger.log(LogLevel::WARNING, "Strong Pruning: Found %zu components.\n", final_components_.size());

            for (size_t comp_idx = 0; comp_idx < final_components_.size(); ++comp_idx) {
                if(final_components_[comp_idx].empty()) { context.logger.log(LogLevel::INFO, "  Skipping empty component %zu.\n", comp_idx); continue; }
                context.logger.log(LogLevel::WARNING,"Strong Pruning: Processing component %zu (size %zu).\n", comp_idx, final_components_[comp_idx].size());

                PCSTFast::IndexType root_node;
                if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                    root_node = context.root;
                    context.logger.log(LogLevel::INFO,"  Using designated root %d for component %zu.\n", root_node, comp_idx);
                } else {
                    root_node = find_best_component_root(context, static_cast<PCSTFast::IndexType>(comp_idx));
                    context.logger.log(LogLevel::INFO,"  Using best root %d for component %zu.\n", root_node, comp_idx);
                }

                if(root_node != PCSTFast::kInvalidIndex && static_cast<size_t>(root_node) < context.prizes.size()) {
                    context.logger.log(LogLevel::INFO,"  Running final strong_pruning_dfs from %d (MarkDeleted=true).\n", root_node);

                     if(static_cast<PCSTFast::IndexType>(comp_idx) == root_component_index_) {
                            context.logger.log(LogLevel::DEBUG,"  Resetting strong pruning parent/payoff for root component %zu nodes before final DFS.\n", comp_idx);
                         for(PCSTFast::IndexType n : final_components_[comp_idx]) {
                                if(static_cast<size_t>(n) < strong_pruning_parent_.size()) strong_pruning_parent_[n] = {PCSTFast::kInvalidIndex, 0.0};
                                if(static_cast<size_t>(n) < strong_pruning_payoff_.size()) strong_pruning_payoff_[n] = 0.0;
                         }
                     }

                    strong_pruning_dfs(context, root_node, true);
                } else {
                    context.logger.log(LogLevel::WARNING,"Warning: Skipping strong pruning for component %zu due to invalid root (%d).\n", comp_idx, root_node);
                }
            }

            phase3_result_local_.clear();
            phase3_result_local_.reserve(phase2_result_local_.size());
            context.logger.log(LogLevel::WARNING, "Strong pruning component processing complete. Filtering phase 2 edges based on deleted nodes...\n");
            for(PCSTFast::IndexType edge_idx : phase2_result_local_) {
                 if (static_cast<size_t>(edge_idx) >= context.edges.size()) continue;
                const auto& edge = context.edges[edge_idx];
                 bool u_del = (static_cast<size_t>(edge.first) >= node_deleted_.size() || node_deleted_[edge.first]);
                 bool v_del = (static_cast<size_t>(edge.second) >= node_deleted_.size() || node_deleted_[edge.second]);

                if (!u_del && !v_del) {
                    phase3_result_local_.push_back(edge_idx);
                } else {
                     context.logger.log(LogLevel::DEBUG,"  Strong: Edge %d (%d,%d) removed due to deleted endpoint(s).\n", edge_idx, edge.first, edge.second);
                }
            }
            context.logger.log(LogLevel::INFO, "StrongPruner::run_strong_pruning Exit. Final edge count %zu.\n", phase3_result_local_.size());
        }

        void StrongPruner::prune(const PruningContext& context,
                std::vector<PCSTFast::IndexType>& result_nodes,
                std::vector<PCSTFast::IndexType>& result_edges) {
            context.logger.log(LogLevel::INFO, "Pruning: Strong. Setting up...\n");
            setup(context);
            context.logger.log(LogLevel::INFO, "Pruning: Running Strong pruning logic...\n");
            run_strong_pruning(context);
            result_edges = phase3_result_local_;
            context.logger.log(LogLevel::WARNING, "Strong pruning complete. Building final node set...\n");
            build_pruned_node_set(context, result_nodes);
            context.logger.log(LogLevel::INFO, "Final Result (Strong Pruning): Nodes=%zu, Edges=%zu\n", result_nodes.size(), result_edges.size());
        }

    }
}
