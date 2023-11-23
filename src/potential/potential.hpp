#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include "potential/params.hpp"
#include "potential/params_base.hpp"
#include "solver/channel.hpp"
#include "util/types.hpp"

#include <complex>
#include <memory>
#include <vector>

namespace osiris {
/// @brief CRTP interface for a local potential.
/// Local implied a potential V is diagonal in coordinate space,
/// e.g. <r|V|rp> = V(r) delta(r-rp)
template <class T> struct Potential {
  using CH = Channel;
  using AM = Channel::FermionSpinOrbitCoupling;
  using E = Channel::Energetics;

  static constexpr bool is_local = true;
  static constexpr bool is_symmetric = true;

  cmpl eval(real r, const CH& ch, const E& e, const AM& am) const {
    return (static_cast<const T*>(this))->eval_impl(r, ch, e, am);
  }
  cmpl eval_reduced(real r, const CH& ch, const E& e, const AM& am) const {
    return r * eval(r, ch, e, am);
  }
  cmpl operator()(real r, const CH& ch, const E& e, const AM& am) const {
    return eval(r, ch, e, am);
  }
};

/// @brief CRTP interface for a non-local potential that is symmetric in the
/// arguments e.g. V(r1,r2) = V(r2,r1)
template <class T> struct NonlocalPotential {
  using CH = Channel;
  using AM = Channel::FermionSpinOrbitCoupling;
  using E = Channel::Energetics;
  static constexpr bool is_local = false;
  static constexpr bool is_symmetric = true;

  cmpl eval(real r, real rp, const CH& ch, const E& e, const AM& am) const {
    return (static_cast<const T*>(this))->eval_impl(r, rp, ch, e, am);
  }
  cmpl
  eval_reduced(real r, real rp, const CH& ch, const E& e, const AM& am) const {
    return eval(r, rp, ch, e, am) * r * rp;
  }
  cmpl eval_reduced(real r, const CH& ch, const E& e, const AM& am) const {
    return r * r * eval(r, r, ch, e, am);
  }
  cmpl eval(real r, const CH& ch, const E& e, const AM& am) const {
    return eval(r, r, ch, e, am);
  }
  cmpl
  operator()(real r, real rp, const CH& ch, const E& e, const AM& am) const {
    return eval(r, rp, ch, e, am);
  }
};

/// @brief Common phenomenological potential form used for central potentials
class WoodsSaxon : public Potential<WoodsSaxon> {
private:
  cmpl V;
  real R, a;

public:
  WoodsSaxon(real V, real R, real a) : V(V), R(R), a(a){};
  WoodsSaxon(cmpl V, real R, real a) : V(V), R(R), a(a){};

  WoodsSaxon(const WoodsSaxon& rhs) = default;
  WoodsSaxon(WoodsSaxon&& rhs) = default;
  WoodsSaxon() = default;

  cmpl eval_impl(real r, const CH&, const E&, const AM&) const {
    return V / (1. + exp((r - R) / a));
  };
};

/// @brief Common phenomenological potential form used for surface peakes
/// potentials. The derivative in r of a Woods-Saxon
class DerivWoodsSaxon : public Potential<DerivWoodsSaxon> {
private:
  cmpl V;
  real R, a;

public:
  DerivWoodsSaxon(real V, real R, real a) : V(V), R(R), a(a){};
  DerivWoodsSaxon(cmpl V, real R, real a) : V(V), R(R), a(a){};

  DerivWoodsSaxon(const DerivWoodsSaxon& rhs) = default;
  DerivWoodsSaxon(DerivWoodsSaxon&& rhs) = default;
  DerivWoodsSaxon() = default;

  cmpl eval_impl(real r, const CH&, const E&, const AM&) const {
    const auto y = exp((r - R)) / a;
    return -V / a * (y / ((1. + y) * (1. + y)));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
class Thomas : public Potential<Thomas> {
private:
  cmpl V;
  real R, a;

public:
  Thomas(real V, real R, real a) : V(V), R(R), a(a){};
  Thomas(cmpl V, real R, real a) : V(V), R(R), a(a){};

  Thomas(const Thomas& rhs) = default;
  Thomas(Thomas&& rhs) = default;
  Thomas() = default;

  cmpl eval_impl(real r, const CH& ch, const E& e, const AM& am) const {
    return DerivWoodsSaxon(V, R, a).eval(r, ch, e, am) / r;
  };
};

class Gaussian : public Potential<Gaussian> {
private:
  cmpl V;
  real R, sigma;

public:
  Gaussian(real V, real R, real sigma) : V(V), R(R), sigma(sigma){};
  Gaussian(cmpl V, real R, real sigma) : V(V), R(R), sigma(sigma){};

  Gaussian(const Gaussian& rhs) = default;
  Gaussian(Gaussian&& rhs) = default;
  Gaussian() = default;

  cmpl eval_impl(real r, const CH&, const E&, const AM&) const {
    return exp((r - R) * (r - R) / (sigma * sigma));
  };
};

template <size_t N> class NGaussian : public Potential<NGaussian<N>> {
private:
  std::array<Gaussian, N> gaussians;

public:
  using typename Potential<NGaussian<N>>::CH;
  using typename Potential<NGaussian<N>>::E;
  using typename Potential<NGaussian<N>>::AM;

  NGaussian(std::array<Gaussian, N>&& gaussians)
      : gaussians(std::move(gaussians)){};

  NGaussian(const NGaussian& rhs) = default;
  NGaussian(NGaussian&& rhs) = default;

  cmpl eval_impl(real r, const CH& ch, const E& e, const AM& am) const {
    cmpl v = 0;
    for (unsigned int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r, ch, e, am);
    }
    return v;
  };
};

using DoubleGaussian = NGaussian<2>;

///@brief Thompson, D. R., LeMere, M., & Tang, Y. C. (1977).
/// Systematic investigation of scattering problems with the resonating-group
/// method. Nuclear Physics A, 286(1), 53-66.
class Minnesota : public DoubleGaussian {
public:
  static Minnesota build_1S0() { return Minnesota(200, -98.15, 1.487, 0.465); }
  static Minnesota build_3S1() { return Minnesota(200, -178, 1.487, 0.639); }

  Minnesota(const Minnesota& rhs) = default;
  Minnesota(Minnesota&& rhs) = default;

  Minnesota(real V01, real V02, real k1, real k2)
      : DoubleGaussian(
            {Gaussian(V01, 0., 1 / sqrt(k1)),
             Gaussian(V02, 0., 1 / sqrt(k2))}) {}
};

class Yukawa : public Potential<Yukawa> {
private:
  real mass_coupling;
  real force_coupling;

public:
  cmpl eval_impl(real r, const CH&, const E&, const AM&) const {
    return -force_coupling * exp(-r / mass_coupling) / r;
  };
  Yukawa(real mass_coupling, real force_coupling)
      : mass_coupling(mass_coupling), force_coupling(force_coupling){};
  Yukawa(real alpha, real m, real g)
      : mass_coupling(alpha * m), force_coupling(g * g){};
};

class SphereWell : public Potential<SphereWell> {
private:
  real R;
  cmpl depth;

public:
  cmpl eval_impl(real r, const CH&, const E&, const AM&) const {
    if (r < R)
      return -depth;
    return 0;
  };

  SphereWell(real R, cmpl depth) : R(R), depth(depth){};
};

/// @brief Common global phenomenological OMP form for a given isotope, A,Z.
template <class Params> class OMP : public Potential<OMP<Params>> {
public:
  struct Terms {

    WoodsSaxon V, W;
    DerivWoodsSaxon Vs, Ws;
    Thomas Vso, Wso;

    Terms(int A, int Z, real erg_lab, const Params& p)
        : V(WoodsSaxon{
              cmpl{p.real_cent_V(Z, A, erg_lab), 0},
              p.real_cent_r(Z, A, erg_lab), p.real_cent_a(Z, A, erg_lab)}),
          W(WoodsSaxon{
              cmpl{0, p.cmpl_cent_V(Z, A, erg_lab)},
              p.cmpl_cent_r(Z, A, erg_lab), p.cmpl_cent_a(Z, A, erg_lab)}),

          Vs(DerivWoodsSaxon{
              cmpl{p.real_surf_V(Z, A, erg_lab), 0},
              p.real_surf_r(Z, A, erg_lab), p.real_surf_a(Z, A, erg_lab)}),
          Ws(DerivWoodsSaxon{
              cmpl{0, p.cmpl_surf_V(Z, A, erg_lab)},
              p.cmpl_surf_r(Z, A, erg_lab), p.cmpl_surf_a(Z, A, erg_lab)}),

          Vso(Thomas{
              cmpl{p.real_spin_V(Z, A, erg_lab), 0},
              p.real_spin_r(Z, A, erg_lab), p.real_spin_a(Z, A, erg_lab)}),
          Wso(Thomas{
              cmpl{0, p.cmpl_spin_V(Z, A, erg_lab)},
              p.cmpl_spin_r(Z, A, erg_lab), p.cmpl_spin_a(Z, A, erg_lab)}) {}
  };

private:
  /// @brief target
  int A, Z;

  /// @brief Determines the parameterization used for calculating
  /// term depths, radii, and diffusivities
  Params params;

public:
  using typename Potential<OMP<Params>>::CH;
  using typename Potential<OMP<Params>>::E;
  using typename Potential<OMP<Params>>::AM;

  OMP(Isotope isotope) : A(isotope.A), Z(isotope.Z), params(Params()){};
  OMP(Isotope isotope, Params params)
      : A(isotope.A), Z(isotope.Z), params(params){};
  OMP(Isotope isotope, Params&& params)
      : A(isotope.A), Z(isotope.Z), params(params){};

  OMP(const OMP<Params>& rhs) = default;
  OMP(OMP<Params>&& rhs) = default;

  void reset_target(int An, int Zn) {
    A = An;
    Z = Zn;
  }

  cmpl eval_set_terms(
      real r, const Terms& terms, const CH& ch, const E& e, const AM& am) {
    const auto [V, W, Vs, Ws, Vso, Wso] = terms;
    const auto so = am.spin_orbit();

    return V.eval(r, ch, e, am) + Vs.eval(r, ch, e, am) +
           Vso.eval(r, ch, e, am) * so // R
           + W.eval(r, ch, e, am) + Ws.eval(r, ch, e, am) +
           Wso.eval(r, ch, e, am) * so; // Im
  }

  cmpl eval_impl(real r, const CH& ch, const E& e, const AM& am) const {
    const auto terms = Terms(A, Z, e.erg_lab, params);
    return eval_set_terms(r, terms, ch, e, am);
  }
};

/// @brief An arbitrary local potential in r smeared into the off-diagonal
/// by a Gaussian factor in (r-rp), from:
/// Perey, F., and B. Buck.
/// "A non-local potential model for the scattering of neutrons by nuclei."
/// Nuclear Physics 32 (1962): 353-380.
template <class LocalPotential>
class PereyBuck : public NonlocalPotential<PereyBuck<LocalPotential>> {
private:
  LocalPotential local_potential;
  Gaussian non_local_factor;

public:
  using typename Potential<PereyBuck<LocalPotential>>::CH;
  using typename Potential<PereyBuck<LocalPotential>>::E;
  using typename Potential<PereyBuck<LocalPotential>>::AM;

  PereyBuck(LocalPotential potential, real beta)
      : local_potential(std::move(potential)),
        non_local_factor(
            1 / (pow(constants::pi, 3. / 2.)) * beta * beta * beta, 0, beta) {}

  cmpl
  eval_impl(real r, real rp, const CH& ch, const E& e, const AM& am) const {
    return local_potential->eval(r, ch, e, am) *
           non_local_factor.eval((r - rp), ch, e, am);
  };
};

/// @brief Yamaguchi, Yoshio.
/// "Two-nucleon problem when the potential is nonlocal but separable. I."
/// Physical Review 95.6 (1954): 1628.
class Yamaguchi : public NonlocalPotential<Yamaguchi> {
private:
  real alpha;
  real beta;
  real f;

public:
  /// @param mu reduced mass [amu]
  /// @param alpha [fm]^-1
  /// @param beta [fm]^-1
  Yamaguchi(real mu, real alpha, real beta)
      : alpha(alpha), beta(beta),
        f(constants::hbar * constants::hbar * constants::c * constants::c /
          (mu * constants::MeV_per_amu)){};

  ///@brief parameters chosen to reproduce bound state of deuteron and
  /// and neutron-proton triplet scattering length
  Yamaguchi() : alpha(0.2316053), beta(1.3918324), f(41.472){};

  cmpl eval_impl(real r, real rp, const CH&, const E&, const AM&) const {
    return f * 2 * beta * (alpha + beta) * (alpha + beta) *
           exp(-beta * (r + rp));
  };

  real analytic_swave_kmatrix(real k) const {
    const auto a = alpha;
    const auto b = beta;
    const auto d = 2 * (a + b) * (a + b);
    real cot_delta = (a * b * (a + 2 * b) / d +
                      k * k * (a * a + 2 * a * b + 3 * b * b) / (b * d) +
                      k * k * k * k / (b * d)) /
                     k;
    return 1. / cot_delta;
  }
};

} // namespace osiris

#endif
