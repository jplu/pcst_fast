#include "pcst_fast/logger.h"
#include "test_helpers.h"

#include <sstream>
#include <string>

#include "gtest/gtest.h"

using namespace cluster_approx;

class MockLogger : public Logger {
  public:
    std::stringstream log_stream;

  protected:
    void log_impl(LogLevel level, const std::string& message) override {

        log_stream << static_cast<int>(level) << ": " << message << std::endl;
    }
};

TEST(LoggerTest, LevelSetting) {
    StderrLogger logger;

    logger.set_level(LogLevel::INFO);
    EXPECT_EQ(logger.get_level(), LogLevel::INFO);

    logger.set_level(LogLevel::DEBUG);
    EXPECT_EQ(logger.get_level(), LogLevel::DEBUG);

    logger.set_level(LogLevel::FATAL);
    EXPECT_EQ(logger.get_level(), LogLevel::FATAL);
}

TEST(LoggerTest, NullLoggerWorks) {
    test_utils::NullLogger logger;

    EXPECT_NO_THROW(logger.log(LogLevel::TRACE, "Trace message {}", 1));
    EXPECT_NO_THROW(logger.log(LogLevel::INFO, "Info message"));
    EXPECT_NO_THROW(logger.log(LogLevel::ERROR, "Error message"));
}

TEST(LoggerTest, LevelFiltering) {
    MockLogger logger;

    logger.set_level(LogLevel::INFO);

    logger.log(LogLevel::ERROR, "Error should be logged.");
    logger.log(LogLevel::WARNING, "Warning should be logged.");
    logger.log(LogLevel::INFO, "Info should be logged.");
    logger.log(LogLevel::DEBUG, "Debug should NOT be logged.");
    logger.log(LogLevel::TRACE, "Trace should NOT be logged.");

    std::string output = logger.log_stream.str();
    EXPECT_NE(output.find("1: Error should be logged."), std::string::npos);
    EXPECT_NE(output.find("2: Warning should be logged."), std::string::npos);
    EXPECT_NE(output.find("3: Info should be logged."), std::string::npos);
    EXPECT_EQ(output.find("4: Debug should NOT be logged."), std::string::npos);
    EXPECT_EQ(output.find("5: Trace should NOT be logged."), std::string::npos);

    logger.set_level(LogLevel::DEBUG);
    logger.log(LogLevel::DEBUG, "Debug should NOW be logged.");

    output = logger.log_stream.str();
    EXPECT_NE(output.find("4: Debug should NOW be logged."), std::string::npos);

}

TEST(LoggerTest, Formatting) {
    MockLogger logger;
    logger.set_level(LogLevel::INFO);

    int i = 10;
    double d = 3.14;
    std::string s = "test";
    logger.log(LogLevel::INFO, "Log with args: int={}, double={:.3f}, string='{}'", i, d, s);

    std::string output = logger.log_stream.str();

    EXPECT_NE(output.find("Log with args: int=10"), std::string::npos);
    EXPECT_NE(output.find("double=3.140"), std::string::npos);
    EXPECT_NE(output.find("string='test'"), std::string::npos);
}

