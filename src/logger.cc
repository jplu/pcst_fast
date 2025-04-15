#include "pcst_fast/logger.h"

#include <cstdio>
#include <chrono>
#include <iomanip>
#include <sstream>

namespace cluster_approx {

StderrLogger::StderrLogger(LogLevel initial_level) {
    current_level_ = initial_level;
}

const char* StderrLogger::level_to_string(LogLevel level) {
    switch (level) {
    case LogLevel::FATAL:
        return "FATAL";
    case LogLevel::ERROR:
        return "ERROR";
    case LogLevel::WARNING:
        return "WARN ";
    case LogLevel::INFO:
        return "INFO ";
    case LogLevel::DEBUG:
        return "DEBUG";
    case LogLevel::TRACE:
        return "TRACE";
    default:
        return "?????";
    }
}

void StderrLogger::log_impl(LogLevel level, const std::string& message) {
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::ostringstream time_stream;
#ifdef _MSC_VER
    std::tm now_tm;
    localtime_s(&now_tm, &now_c);
    time_stream << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
#else
    std::tm now_tm;
    localtime_r(&now_c, &now_tm);
    time_stream << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
#endif
    time_stream << '.' << std::setfill('0') << std::setw(3) << now_ms.count();

    fprintf(stderr, "[%s] [%s] %s\n",
            time_stream.str().c_str(),
            level_to_string(level),
            message.c_str());
    fflush(stderr);
}

}