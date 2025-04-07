#include "pruner_factory.h"

#include "ipruner.h"
#include "no_pruner.h"
#include "simple_pruner.h"
#include "gw_pruner.h"
#include "strong_pruner.h"
#include "connect_final_pruner.h"

#include <memory>
#include <stdexcept>

namespace cluster_approx {
    namespace internal {

        std::unique_ptr<IPruner> create_pruner(PCSTFast::PruningMethod method) {
            switch (method) {
                case PCSTFast::PruningMethod::kNoPruning:
                    return std::make_unique<NoPruner>();
                case PCSTFast::PruningMethod::kSimplePruning:
                    return std::make_unique<SimplePruner>();
                case PCSTFast::PruningMethod::kGWPruning:
                    return std::make_unique<GWPruner>();
                case PCSTFast::PruningMethod::kStrongPruning:
                    return std::make_unique<StrongPruner>();
                case PCSTFast::PruningMethod::kConnectFinalComponents:
                    return std::make_unique<ConnectFinalPruner>();
                case PCSTFast::PruningMethod::kUnknownPruning:
                default:
                    throw std::invalid_argument("Unsupported or unknown pruning method provided to factory.");
            }
        }

    }
}
