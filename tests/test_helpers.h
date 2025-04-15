#pragma once // Use pragma once

#include <vector>
#include <cstdio>
#include <algorithm>

#include "gtest/gtest.h"

template <typename T, size_t N>
[[nodiscard]] constexpr T* begin(T(&arr)[N]) noexcept {
    return &arr[0];
}
template <typename T, size_t N>
[[nodiscard]] constexpr T* end(T(&arr)[N]) noexcept {
    return &arr[0] + N;
}


/**
 * @brief Compares two vectors element by element after sorting them.
 * Asserts that sizes are equal and elements match. Used for comparing results like node/edge sets.
 * @tparam T The type of elements in the vectors. Must support operator== and operator<.
 * @param expected_result The vector containing the expected elements.
 * @param result The vector containing the actual elements obtained.
 */
template <typename T>
void CheckResult(const std::vector<T>& expected_result, const std::vector<T>& result) {
    std::vector<T> sorted_expected = expected_result;
    std::vector<T> sorted_result = result;

    std::sort(sorted_expected.begin(), sorted_expected.end());
    std::sort(sorted_result.begin(), sorted_result.end());

    ASSERT_EQ(sorted_expected.size(), sorted_result.size()) << "Result vectors have different sizes.";

    for (size_t i = 0; i < sorted_expected.size(); ++i) {
        EXPECT_EQ(sorted_expected[i], sorted_result[i])
                << "Vectors differ at index " << i << " after sorting.";
    }
}

#include "pcst_fast/logger.h"

namespace cluster_approx {
namespace test_utils {

/**
 * @brief A logger implementation that does nothing. Useful for tests.
 */
class NullLogger : public Logger {
  public:
    NullLogger() {
        set_level(LogLevel::TRACE);
    }
  protected:
    void log_impl(LogLevel, const std::string&) override {
    }
};

}
}