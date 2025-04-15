#pragma once

#include <string_view>
#include <string>
#include <stdexcept>
#include <format>

namespace cluster_approx {

/**
 * @brief Defines the severity levels for logging messages.
 */
enum class LogLevel {
    FATAL = 0,
    ERROR,
    WARNING,
    INFO,
    DEBUG,
    TRACE
};

/**
 * @brief Abstract interface for logging messages.
 *
 * Implementations of this interface handle the actual output of log messages.
 * The logger can be configured with a minimum level to filter messages.
 */
class Logger {
  public:
    virtual ~Logger() = default;

    /**
     * @brief Logs a message if its level is greater than or equal to the configured minimum level.
     * @param level The severity level of the message.
     * @param fmt A std::format string.
     * @param args Arguments for the format string.
    */
    template<typename... Args>
    void log(LogLevel level, std::format_string<Args...> fmt, Args&&... args) {
        if (level <= current_level_) {
            try {
                log_impl(level, std::format(fmt, std::forward<Args>(args)...));
            } catch (const std::format_error& fe) {
                log_impl(LogLevel::ERROR, std::string("Logging format error: ") + fe.what());
            } catch (...) {
                log_impl(LogLevel::ERROR, "Unknown error during logging formatting.");
            }
        }
    }

    /**
     * @brief Sets the minimum logging level. Messages below this level will be ignored.
     * @param level The minimum LogLevel to output.
     */
    void set_level(LogLevel level) {
        current_level_ = level;
    }

    /**
     * @brief Gets the current minimum logging level.
     * @return The current LogLevel.
     */
    [[nodiscard]] LogLevel get_level() const {
        return current_level_;
    }

  protected:
    /**
     * @brief Abstract method to be implemented by subclasses for actual log output.
     * @param level The severity level of the message.
     * @param message The fully formatted message string.
     */
    virtual void log_impl(LogLevel level, const std::string& message) = 0;

    LogLevel current_level_ = LogLevel::INFO;
};

/**
 * @brief A simple logger implementation that writes to stderr.
 */
class StderrLogger : public Logger {
  public:
    /**
     * @brief Constructs a StderrLogger.
     * @param initial_level The initial minimum logging level. Defaults to INFO.
     */
    explicit StderrLogger(LogLevel initial_level = LogLevel::INFO);

  protected:
    void log_impl(LogLevel level, const std::string& message) override;

  private:
    static const char* level_to_string(LogLevel level);
};

}