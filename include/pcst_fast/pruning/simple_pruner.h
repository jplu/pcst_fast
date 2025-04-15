#pragma once

#include "pcst_fast/pcst_interfaces.h"

namespace cluster_approx {
namespace pruning {

/**
 * @brief Implements the "simple pruning" strategy.
 * Returns the intermediate result from the core algorithm directly.
 * This is functionally equivalent to NoPruner given the definition of CoreAlgorithmResult.
 */
class SimplePruner final : public IPruner {
  public:
    /**
     * @brief Applies the "simple pruning" strategy.
     * @param input Data structure containing the graph and the result from the core algorithm.
     * @return A PruningResult containing the intermediate edges and corresponding nodes.
     * @throws std::exception or derived classes on error (e.g., from utility functions).
     */
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;
};

}
}