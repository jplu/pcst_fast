#include "pcst_fast/pcst_types.h"
#include <string>
#include <algorithm>
#include <cctype>

namespace cluster_approx {

PruningMethod parse_pruning_method(const std::string& input) {
    std::string input_lower;
    input_lower.resize(input.size());
    std::transform(input.begin(), input.end(), input_lower.begin(), ::tolower);

    if (input_lower == "none") {
        return PruningMethod::kNone;
    } else if (input_lower == "simple") {
        return PruningMethod::kSimple;
    } else if (input_lower == "gw") {
        return PruningMethod::kGW;
    } else if (input_lower == "strong") {
        return PruningMethod::kStrong;
    } else {
        return PruningMethod::kUnknown;
    }
}

}