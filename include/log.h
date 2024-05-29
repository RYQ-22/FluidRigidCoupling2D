#ifndef LOG_H
#define LOG_H

#include "config.h"

namespace backend {

void Assert(const bool condition, const std::string& location, const std::string& message);

}

#endif