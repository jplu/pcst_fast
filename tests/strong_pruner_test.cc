#include "pcst_fast/pruning/strong_pruner.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "test_helpers.h"

#include <vector>
#include <utility>
#include <span>

#include "gtest/gtest.h"

using namespace cluster_approx;

class StrongPrunerTest : public ::testing::Test {
  protected:
    test_utils::NullLogger logger;

    std::vector<std::pair<NodeId, NodeId>> edges_vec = {{0, 1}, {1, 2}, {2, 3}};
    std::vector<double> prizes_vec = {10.0, 10.0, 10.0, 10.0};
    std::vector<double> costs_vec = {1.0, 1.0, 1.0};
    GraphData graph;
    CoreAlgorithmResult core_result;

    StrongPrunerTest() :
        graph{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId}
    {}

    void SetUp() override {
        core_result.initial_node_filter.assign(prizes_vec.size(), true);
        core_result.phase1_edges = {0, 1, 2};
        core_result.edge_inactive_merge_event_ids.assign(3, kInvalidEventId);
        core_result.inactive_merge_events.clear();
    }
};

TEST_F(StrongPrunerTest, NoPruningNeeded) {
    // High prizes, low costs, everything kept
    prizes_vec = {10.0, 10.0, 10.0, 10.0};
    costs_vec = {1.0, 1.0, 1.0};
    graph.prizes = std::span(prizes_vec);
    graph.costs = std::span(costs_vec);

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0, 1, 2};
    std::vector<NodeId> expected_nodes = {0, 1, 2, 3};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}

TEST_F(StrongPrunerTest, PruneTerminalEdge) {
    // Last node (3) has low prize (1.0), edge cost is high (5.0). Net < 0.
    prizes_vec = {10.0, 10.0, 10.0, 1.0};
    costs_vec = {1.0, 1.0, 5.0};
    graph.prizes = std::span(prizes_vec);
    graph.costs = std::span(costs_vec);

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0, 1};
    std::vector<NodeId> expected_nodes = {0, 1, 2};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}

TEST_F(StrongPrunerTest, PruneMiddleEdge) {
    // 0(10) --1--> 1(1) --15--> 2(1) --1--> 3(10)
    // Both 1 and 2 are low prize. The bridge cost 15 is too high for the subtrees.
    prizes_vec = {10.0, 1.0, 1.0, 10.0};
    costs_vec = {1.0, 15.0, 1.0};
    graph.prizes = std::span(prizes_vec);
    graph.costs = std::span(costs_vec);

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {};
    std::vector<NodeId> expected_nodes = {0};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}

TEST_F(StrongPrunerTest, TwoComponents) {
    // Disconnected phase 1 edges: {0,1} and {2,3}
    edges_vec = {{0, 1}, {2, 3}};
    prizes_vec = {10.0, 10.0, 5.0, 5.0};
    costs_vec = {1.0, 1.0};

    graph = GraphData{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId};

    SetUp();
    core_result.initial_node_filter = {true, true, true, true};
    core_result.phase1_edges = {0, 1};

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0, 1};
    std::vector<NodeId> expected_nodes = {0, 1, 2, 3};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}

TEST_F(StrongPrunerTest, TwoComponentsOnePruned) {
    // Disconnected: {0,1} high value, {2,3} low value high cost
    edges_vec = {{0, 1}, {2, 3}};
    prizes_vec = {10.0, 10.0, 1.0, 1.0};
    costs_vec = {1.0, 5.0}; // Edge (2,3) costs 5, total prize 2. Pruned.

    graph = GraphData{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId};

    SetUp();
    core_result.initial_node_filter = {true, true, true, true};
    core_result.phase1_edges = {0, 1};

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0};
    
    // Result logic:
    // {0, 1}: Prize 20, Cost 1. Kept.
    // {2, 3}: Prize 2, Cost 5. 
    //   Subtree {2}: Prize 1, Cost 0. Profit 1.
    //   Subtree {3}: Prize 1, Cost 0. Profit 1.
    //   Subtree {2,3}: Prize 2, Cost 5. Profit -3.
    //   The pruner maximizes payoff, so it keeps node 2 (as an isolated component) 
    //   and prunes the edge to 3. Node 3 is removed.
    //   Nodes surviving: {0, 1, 2}.
    std::vector<NodeId> strict_expected_nodes = {0, 1, 2};
    
    CheckResult(expected_edges, result.edges);
    CheckResult(strict_expected_nodes, result.nodes);
}