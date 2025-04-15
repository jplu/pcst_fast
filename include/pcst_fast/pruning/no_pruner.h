#pragma once

#include "pcst_fast/pcst_interfaces.h"

namespace cluster_approx {
namespace pruning {

/**
 * @brief Implements the "no pruning" strategy.
 * Returns the intermediate result from the core algorithm directly.
 * The nodes are derived from the edges included in the intermediate result
 * and any initially good nodes that might be isolated.
 */
class NoPruner final : public IPruner {
  public:
    /**
     * @brief Applies the "no pruning" strategy.
     * @param input Data structure containing the graph and the result from the core algorithm.
     * @return A PruningResult containing the intermediate edges and corresponding nodes.
     * @throws std::exception or derived classes on error (e.g., from utility functions).
     */
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;
};

}
}