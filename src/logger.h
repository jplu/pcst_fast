#pragma once

#include <string>
#include <functional>
#include <cstdarg>
#include <cstdio>
#include <memory>
#include <mutex>
#include <atomic>

namespace cluster_approx {
    namespace internal {

        enum class LogLevel {
            NONE = -1,
            FATAL = 0,
            ERROR = 1,
            WARNING = 2,
            INFO = 3,
            DEBUG = 4,
            TRACE = 5
        };

        using LogSink = std::function<void(LogLevel level, const std::string& message)>;

        void DefaultLogSink(LogLevel level, const std::string& message);

        class Logger {
            private:
                LogSink sink_;
                std::atomic<LogLevel> current_level_;
                std::mutex sink_mutex_;

                void vlog(LogLevel level, const char* format, va_list args) noexcept;

            public:
                explicit Logger(LogSink sink = DefaultLogSink, LogLevel level = LogLevel::INFO);

                Logger(const Logger&) = delete;
                Logger& operator=(const Logger&) = delete;
                Logger(Logger&&) = delete;
                Logger& operator=(Logger&&) = delete;

                void set_level(LogLevel level) noexcept;

                [[nodiscard]] LogLevel get_level() const noexcept;

                [[nodiscard]] bool is_enabled(LogLevel level) const noexcept;

                void log(LogLevel level, const char* format, ...) const noexcept
                    #if defined(__GNUC__) || defined(__clang__)
                        __attribute__ ((format (printf, 3, 4)))
                    #endif
                    ;

                void log(LogLevel level, const char* format, ...) noexcept
                    #if defined(__GNUC__) || defined(__clang__)
                        __attribute__ ((format (printf, 3, 4)))
                    #endif
                    ;
        };

        inline LogLevel TranslateLegacyLogLevel(int old_level) {
            switch (old_level) {
                case -1: return LogLevel::NONE;
                case 0:  return LogLevel::FATAL;
                case 1:  return LogLevel::ERROR;
                case 2:  return LogLevel::WARNING;
                case 3:  return LogLevel::INFO;
                case 4:  return LogLevel::DEBUG;
                case 5:  return LogLevel::TRACE;
                default:
                    if (old_level < -1) return LogLevel::NONE;
                    else return LogLevel::TRACE;
            }
        }
    }
}
