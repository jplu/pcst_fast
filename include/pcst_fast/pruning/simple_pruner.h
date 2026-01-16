#pragma once

#include "pcst_fast/pcst_interfaces.h"

namespace cluster_approx {
namespace pruning {

class SimplePruner final : public IPruner {
  public:
    [[nodiscard]] PruningResult prune(const PruningInput& input) override;
};

}
}