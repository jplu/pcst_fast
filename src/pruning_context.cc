#include "pruning_context.h"

#include <cstdio>
#include <string>
#include <cstdarg>
#include <memory>
#include <stdexcept>

namespace cluster_approx {
    namespace internal {

        void PruningContext::log(int level, const char* format, ...) const {
            if (logger && level <= verbosity_level) {
                constexpr int BUFFER_SIZE = 1024;
                char buffer[BUFFER_SIZE];
                va_list args;
                va_start(args, format);
                int needed = vsnprintf(buffer, BUFFER_SIZE, format, args);
                va_end(args);

                if (needed < 0) {
                    logger(level, "[Pruner Logging Error: vsnprintf failed]");
                } else if (static_cast<size_t>(needed) >= BUFFER_SIZE) {
                    buffer[BUFFER_SIZE - 1] = '\0';
                    if (BUFFER_SIZE > 4) { snprintf(buffer + BUFFER_SIZE - 5, 5, "..."); }
                    logger(level, std::string(buffer) + "(TRUNCATED)");
                } else {
                    logger(level, std::string(buffer));
                }
            }
        }

    }
}
