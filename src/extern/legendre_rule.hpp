#include "util/types.hpp"

#include <vector>

namespace gl {

struct zero_crossing {
  osiris::real abscissa;
  osiris::real weight;
};

std::vector<zero_crossing> generate_gauss_legendre_quadrature(int order);

} // namespace gl
