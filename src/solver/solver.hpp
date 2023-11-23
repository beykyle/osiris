#ifndef SOLVER_HEADER
#define SOLVER_HEADER

#include "extern/legendre_rule.hpp"
#include "potential/potential.hpp"
#include "solver/channel.hpp"
#include "util/asymptotics.hpp"
#include "util/config.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

namespace osiris {

class RMatrixKernel {
private:
  const int nbasis{};
  const int nchannels{};
  std::vector<gl::zero_crossing> quadrature;

public:
  RMatrixKernel(int nbasis, int nchannels)
      : nbasis(nbasis), nchannels(nchannels),
        quadrature(gl::generate_gauss_legendre_quadrature(nbasis)) {}
};

} // namespace osiris

#endif
