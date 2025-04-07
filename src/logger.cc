#include "logger.h"
#include <iostream>
#include <stdexcept>

namespace cluster_approx {
    namespace internal {

        void DefaultLogSink(LogLevel level, const std::string& message) {
            switch (level) {
                case LogLevel::FATAL:   std::cerr << "[FATAL] ";   break;
                case LogLevel::ERROR:   std::cerr << "[ERROR] ";   break;
                case LogLevel::WARNING: std::cerr << "[WARNING] "; break;
                case LogLevel::INFO:    std::cerr << "[INFO] ";    break;
                case LogLevel::DEBUG:   std::cerr << "[DEBUG] ";   break;
                case LogLevel::TRACE:   std::cerr << "[TRACE] ";   break;
                case LogLevel::NONE:    return;
            }
            std::cerr << message << std::flush;
        }

        Logger::Logger(LogSink sink, LogLevel level)
            : sink_(std::move(sink)), current_level_(level) {
            if (!sink_) {
                sink_ = DefaultLogSink;
            }
        }

        void Logger::set_level(LogLevel level) noexcept {
            current_level_.store(level, std::memory_order_release);
        }

        LogLevel Logger::get_level() const noexcept {
            return current_level_.load(std::memory_order_acquire);
        }

        bool Logger::is_enabled(LogLevel level) const noexcept {
            return level <= current_level_.load(std::memory_order_acquire) && level != LogLevel::NONE;
        }

        void Logger::vlog(LogLevel level, const char* format, va_list args) noexcept {
            if (!is_enabled(level)) {
                return;
            }

            constexpr int BUFFER_SIZE = 1024;
            char buffer[BUFFER_SIZE];

            va_list args_copy;
            va_copy(args_copy, args);

            int needed = vsnprintf(buffer, BUFFER_SIZE, format, args_copy);

            va_end(args_copy);

            std::string message_str;
            if (needed < 0) {
                message_str = "[Logger Formatting Error]";
            } else if (static_cast<size_t>(needed) >= BUFFER_SIZE) {
                buffer[BUFFER_SIZE - 1] = '\0';
                if (BUFFER_SIZE > 4) {
                    snprintf(buffer + BUFFER_SIZE - 5, 5, "...");
                }
                message_str = std::string(buffer) + "(TRUNCATED)";
            } else {
                message_str = std::string(buffer);
            }

            try {
                sink_(level, message_str);
            } catch (const std::exception& e) {
                try {
                    std::cerr << "[FATAL] Logger sink callback failed: " << e.what() << std::endl;
                } catch (...) {
                }
            } catch (...) {
                try {
                    std::cerr << "[FATAL] Logger sink callback failed with unknown exception." << std::endl;
                } catch (...) { }
            }
        }

        void Logger::log(LogLevel level, const char* format, ...) const noexcept {
            if (!is_enabled(level)) {
                return;
            }

            std::unique_lock<std::mutex> lock(const_cast<Logger*>(this)->sink_mutex_);

            va_list args;
            va_start(args, format);
            const_cast<Logger*>(this)->vlog(level, format, args);
            va_end(args);
        }

        void Logger::log(LogLevel level, const char* format, ...) noexcept {
            if (!is_enabled(level)) {
                return;
            }
            std::unique_lock<std::mutex> lock(sink_mutex_);
            va_list args;
            va_start(args, format);
            vlog(level, format, args);
            va_end(args);
        }

    }
}
