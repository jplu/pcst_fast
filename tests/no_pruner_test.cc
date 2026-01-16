#include "pcst_fast/pruning/no_pruner.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "test_helpers.h"

#include <vector>
#include <utility>
#include <span>

#include "gtest/gtest.h"

using namespace cluster_approx;

TEST(NoPrunerTest, ReturnsIntermediateResult) {
    test_utils::NullLogger logger;
    size_t num_nodes = 5;
    std::vector<std::pair<NodeId, NodeId>> edges_vec = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    std::vector<double> prizes_vec(num_nodes, 1.0);
    std::vector<double> costs_vec = {1.0, 1.0, 1.0, 1.0};

    auto edges_span = std::span<const std::pair<NodeId, NodeId>>(edges_vec);
    auto prizes_span = std::span<const double>(prizes_vec);
    auto costs_span = std::span<const double>(costs_vec);

    GraphData graph{edges_span, prizes_span, costs_span, kInvalidNodeId};

    CoreAlgorithmResult core_result;
    core_result.phase1_edges = {1, 2};
    core_result.initial_node_filter = {false, true, true, true, false};

    PruningInput input{graph, core_result, &logger};

    pruning::NoPruner pruner;
    PruningResult result = pruner.prune(input);

    std::vector<EdgeId> expected_edges = {1, 2};
    CheckResult(expected_edges, result.edges);

    // Nodes 1, 2, 3 are kept because they are in the edges or marked true
    std::vector<NodeId> expected_nodes = {1, 2, 3};
    CheckResult(expected_nodes, result.nodes);
}