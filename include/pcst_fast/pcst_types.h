#pragma once

#include <string>
#include <cstdint>

#if defined(__GNUC__) || defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace cluster_approx {

using NodeId = int32_t;
using EdgeId = int32_t;
using EdgePartId = int32_t;
using ClusterId = int32_t;
using EventId = int32_t;

constexpr NodeId kInvalidNodeId = -1;
constexpr EdgeId kInvalidEdgeId = -1;
constexpr EdgePartId kInvalidEdgePartId = -1;
constexpr ClusterId kInvalidClusterId = -1;
constexpr EventId kInvalidEventId = -1;

enum class PruningMethod {
    kNone = 0,
    kSimple,
    kGW,
    kStrong,
    kUnknown
};

[[nodiscard]] PruningMethod parse_pruning_method(const std::string& input);

}