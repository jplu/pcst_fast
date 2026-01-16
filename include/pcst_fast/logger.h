#pragma once

#include <string_view>
#include <string>
#include <stdexcept>
#include <format>
#include <utility>

namespace cluster_approx {

enum class LogLevel {
    FATAL = 0,
    ERROR,
    WARNING,
    INFO,
    DEBUG,
    TRACE
};

class Logger {
  public:
    virtual ~Logger() = default;

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

    void set_level(LogLevel level) {
        current_level_ = level;
    }

    [[nodiscard]] LogLevel get_level() const {
        return current_level_;
    }

  protected:
    virtual void log_impl(LogLevel level, const std::string& message) = 0;

    LogLevel current_level_ = LogLevel::INFO;
};

class StderrLogger : public Logger {
  public:
    explicit StderrLogger(LogLevel initial_level = LogLevel::INFO);

  protected:
    void log_impl(LogLevel level, const std::string& message) override;

  private:
    static const char* level_to_string(LogLevel level);
};

}