#include "pcst_fast/pruning/gw_pruner.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "pcst_fast/pcst_core_internals.h"
#include "test_helpers.h"

#include <vector>
#include <utility>
#include <span>

#include "gtest/gtest.h"

using namespace cluster_approx;

class GWPrunerTest : public ::testing::Test {
  protected:
    test_utils::NullLogger logger;
    std::vector<std::pair<NodeId, NodeId>> edges_vec = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
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
        core_result.phase1_edges.clear();
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

// Tests for complex active-inactive logic are better handled in EndToEndTest 
// where the 'necessary' flags are set correctly by the core algorithm's structure,
// rather than mocking the fragile internal state here.