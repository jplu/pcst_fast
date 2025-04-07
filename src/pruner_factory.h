#pragma once

#include <memory>
#include "pcst_fast.h"

namespace cluster_approx {
    namespace internal {

        class IPruner;

        std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method);

    }
}
