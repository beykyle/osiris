#ifndef CONFIG_HEADER
#define CONFIG_HEADER

#include "util/types.hpp"

namespace osiris {
/// @brief The maximum number of orbital angular momentum values to consider
/// for a single channel
constexpr int MAXL = 30;
/// @brief Relative error between iterations indicating convergence
constexpr real REL_ERR_EPS = 1.0e-8;
} // namespace osiris

#endif
