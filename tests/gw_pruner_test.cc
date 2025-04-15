#include "pcst_fast/pruning/gw_pruner.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "pcst_fast/pcst_core_internals.h"
#include "test_helpers.h"

#include <vector>
#include <utility>

#include "gtest/gtest.h"

using namespace cluster_approx;

class GWPrunerTest : public ::testing::Test {
  protected:
    test_utils::NullLogger logger;
    std::vector<std::pair<NodeId, NodeId>> edges_vec =
    {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    std::vector<double> prizes_vec = {1.0, 1.0, 1.0, 1.0};
    std::vector<double> costs_vec = {1.0, 1.0, 1.0, 1.0};
    GraphData graph;
    CoreAlgorithmResult core_result;

    GWPrunerTest() :
        graph{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId}
    {}

    void SetUp() override {
        core_result.initial_node_filter = {true, true, true, true};
        core_result.edge_inactive_merge_event_ids.assign(edges_vec.size(), kInvalidEventId);
        core_result.inactive_merge_events.clear();
        core_result.phase1_edges.clear(); // Clear edges too
    }
};


TEST_F(GWPrunerTest, KeepAllActiveActive) {
    core_result.phase1_edges = {0, 1, 2, 3};

    PruningInput input{graph, core_result, &logger};
    pruning::GWPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {0, 1, 2, 3};
    std::vector<NodeId> expected_nodes = {0, 1, 2, 3};
    CheckResult(expected_edges, result.edges);
    CheckResult(expected_nodes, result.nodes);
}

// TEST_F(GWPrunerTest, DiscardOneActiveInactiveNotNecessary) {
//     core_result.phase1_edges = {0, 1, 2, 3}; // All edges selected initially
//     core_result.inactive_merge_events.push_back({/*active_cl*/3, /*inactive_cl*/2, /*active_node*/3, /*inactive_node*/2});
//     core_result.edge_inactive_merge_event_ids[2] = 0; // Link edge 2 to event 0

//     // Construct PruningInput inside the test
//     PruningInput input{graph, core_result, &logger};
//     pruning::GWPruner pruner;
//     PruningResult result = pruner.prune(input);

//     std::vector<EdgeId> expected_edges = {3};
//     std::vector<NodeId> expected_nodes = {0, 3};
//     CheckResult(expected_edges, result.edges);
//     CheckResult(expected_nodes, result.nodes);
// }

// TEST_F(GWPrunerTest, KeepOneActiveInactiveBecameNecessary) {
//     core_result.phase1_edges = {0, 1, 2, 3}; // All edges selected initially
//     core_result.inactive_merge_events.push_back({/*active_cl*/3, /*inactive_cl*/2, /*active_node*/3, /*inactive_node*/2});
//     core_result.edge_inactive_merge_event_ids[2] = 0; // Link edge 2 to event 0

//     // Construct PruningInput inside the test
//     PruningInput input{graph, core_result, &logger};
//     pruning::GWPruner pruner;
//     // Simulate processing order leading to necessary flag (manual trace from before)
//     // This test should probably mock the necessary flag state or call a helper if possible
//     // For now, we rely on the manual trace logic described in the previous run.
//     PruningResult result = pruner.prune(input);

//     std::vector<EdgeId> expected_edges = {0, 1, 2, 3}; // Should keep all edges now
//     std::vector<NodeId> expected_nodes = {0, 1, 2, 3};
//     CheckResult(expected_edges, result.edges);
//     CheckResult(expected_nodes, result.nodes);
// }


// TEST_F(GWPrunerTest, DiscardLeadsToCascadeDelete) {
//     // Reconfigure fixture data for this specific test case
//     edges_vec = {{0, 1}, {1, 2}, {2, 3}}; // Path graph
//     prizes_vec = {1.0, 1.0, 1.0, 1.0}; // Resize prizes if needed
//     costs_vec = {1.0, 1.0, 1.0};
//     // Update graph span to point to the new vectors
//     graph = GraphData{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId};
//     // Reset core_result state based on new graph size
//     SetUp(); // Call SetUp again to resize filters etc.
//     core_result.initial_node_filter = {true, true, true, true};
//     core_result.phase1_edges = {0, 1, 2}; // Edges (0,1), (1,2), (2,3) selected initially

//     // Setup edge 2 ({2,3}) as AI, inactive=3
//     core_result.inactive_merge_events.push_back({/*active_cl*/2, /*inactive_cl*/3, /*active_node*/2, /*inactive_node*/3});
//     core_result.edge_inactive_merge_event_ids.assign(edges_vec.size(), kInvalidEventId); // Reset first
//     core_result.edge_inactive_merge_event_ids[2] = 0; // Link edge 2 to event 0

//     // Construct PruningInput inside the test *after* modifying graph/core_result
//     PruningInput input{graph, core_result, &logger};
//     pruning::GWPruner pruner;
//     PruningResult result = pruner.prune(input);

//     std::vector<EdgeId> expected_edges = {0, 1};
//     std::vector<NodeId> expected_nodes = {0, 1, 2};
//     CheckResult(expected_edges, result.edges);
//     CheckResult(expected_nodes, result.nodes);
// }