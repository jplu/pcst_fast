#include "pcst_fast/pruning/strong_pruner.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "test_helpers.h"

#include <vector>
#include <utility>
#include <cmath>

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
        core_result.edge_inactive_merge_event_ids.clear();
        core_result.inactive_merge_events.clear();
    }
};

TEST_F(StrongPrunerTest, NoPruningNeeded) {
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

    edges_vec = {{0, 1}, {2, 3}};
    prizes_vec = {10.0, 10.0, 1.0, 1.0};
    costs_vec = {1.0, 5.0};

    graph = GraphData{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId};

    SetUp();
    core_result.initial_node_filter = {true, true, true, true};
    core_result.phase1_edges = {0, 1};

    PruningInput input{graph, core_result, &logger};
    pruning::StrongPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0};
    std::vector<NodeId> expected_nodes = {0, 1, 2};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}