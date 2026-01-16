#include "pcst_fast/pcst_core_algorithm.h"
#include "pcst_fast/pcst_interfaces.h"
#include "pcst_fast/pcst_types.h"
#include "pcst_fast/logger.h"
#include "test_helpers.h"

#include <vector>
#include <utility>
#include <span>

#include "gtest/gtest.h"

using namespace cluster_approx;

class PCSTCoreAlgorithmTest : public ::testing::Test {
  protected:
    test_utils::NullLogger logger;
};

TEST_F(PCSTCoreAlgorithmTest, ConstructorThrowsOnInvalidTargetClusters) {
    std::vector<std::pair<NodeId, NodeId>> edges_vec = {{0, 1}};
    std::vector<double> prizes_vec = {1.0, 1.0};
    std::vector<double> costs_vec = {1.0};
    GraphData graph{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), 0};

    // Rooted problem must have target_clusters = 0 (handled internally usually, but constructor enforces 0 input if rooted)
    // Wait, implementation: if (root != -1 && target != 0) throw.
    EXPECT_THROW(PCSTCoreAlgorithm(graph, 1, &logger), std::invalid_argument);

    EXPECT_NO_THROW(PCSTCoreAlgorithm(graph, 0, &logger));
}

TEST_F(PCSTCoreAlgorithmTest, ConstructorThrowsOnNegativePrize) {
    std::vector<std::pair<NodeId, NodeId>> edges_vec = {{0, 1}};
    std::vector<double> prizes_vec = {1.0, -1.0};
    std::vector<double> costs_vec = {1.0};
    GraphData graph{std::span(edges_vec), std::span(prizes_vec), std::span(costs_vec), kInvalidNodeId};

    EXPECT_THROW(PCSTCoreAlgorithm(graph, 1, &logger), std::invalid_argument);
}