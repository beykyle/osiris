#ifndef CHANNEL_HEADER
#define CHANNEL_HEADER

#include "util/asymptotics.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <vector>

namespace osiris {

enum class Parity : bool { odd, even };
enum class Polarization : bool { up, down };

/// @brief  all information relevant to a specific initial/final state
struct Channel {

  /// @brief coupled AM states of the channel. This type is meant to be mutable,
  /// with l and J2 updated according to set_spin
  struct FermionSpinOrbitCoupling {
    /// @brief 2j+1, j being the total angular momentum of channel
    /// (e.g. dimension of the representation of SU(2) denoted by j)
    int j21;
    /// @brief orbital angular momentum [hbar]
    int l;
    /// @brief spatial parity of channel state
    Parity pi;

    Parity update_pi(int l) { return l % 2 == 0 ? Parity::even : Parity::odd; }

    FermionSpinOrbitCoupling(int j21, int l)
        : j21(j21), l(l), pi(update_pi(l)) {
      assert(l >= 0);
      assert(j21 >= 0);
      assert(2 * l + 1 - j21 == -1 or 2 * l + 1 - j21 == 1);
    }

    /// @brief initialize AM to S-Wave state
    FermionSpinOrbitCoupling() : j21(2), l(0), pi(update_pi(l)) {}

    /// @returns projection of spin along axis of orbital ang. mom.; e.g. L * S
    real l_dot_s() const {
      const auto J = static_cast<real>(0.5 * (j21 - 1));
      const auto L = static_cast<real>(l);
      return 0.5 * (J * (J + 1) - L * (L + 1) - 0.5 * (0.5 + 1));
    }

    template <Polarization p>
    static std::vector<FermionSpinOrbitCoupling> generate_couplings(int lmax) {
      std::vector<FermionSpinOrbitCoupling> coupling{};
      int l = 0;
      if constexpr (p == Polarization::up) {
        coupling.resize(lmax);
      } else {
        coupling.resize(lmax - 1);
        l = 1;
      }
      for (; l < lmax; ++l) {
        auto j21 = 2 * l + 2;
        coupling[l] = FermionSpinOrbitCoupling(j21, l);
      }
      return coupling;
    }
  };

  /// @brief things that change w/ energy
  struct Energetics {
    /// @brief CMS energy [MeV]
    real erg_cms;
    /// @brief Bombarding kinetic energy in the lab frame [MeV]
    real erg_lab;
    /// @brief CMS momentum  w/ relativistic correction according to
    /// Eq. 17 of Ingemarsson, A.
    /// "Some notes on optical model calculations at medium energies."
    /// Physica Scripta 9.3 (1974): 156.
    real k;
    /// @brief scattering system reduced mass [MeV]
    real reduced_mass;
    /// @brief hbar^2 * c^2/(2 * reduced mass  * radius) [Mev fm]
    real h2ma;
    /// @brief Sommerfield parameter: sgn(Z_t * Z_p)/(a_B * k) [dimensionless]
    /// ab is the Bohr radius, Z_t and Z_p are the target and product charges
    real sommerfield_param;

    Energetics(real erg_cms, const Channel &ch)
        : erg_cms(erg_cms - ch.threshold),
          erg_lab(erg_cms *
                  (ch.proj_mass * constants::MeV_per_amu +
                   ch.targ_mass * constants::MeV_per_amu) /
                  (ch.targ_mass * constants::MeV_per_amu)),
          k([this, &ch]() -> real {
            using constants::MeV_per_amu;
            using constants::hbar;
            using constants::c;

            // relativistic-corrected wavenumber [1/fm]
            // Eq. 17 of Ingemarsson, 1974
            const auto m1 = ch.proj_mass * MeV_per_amu;
            const auto m2 = ch.targ_mass * MeV_per_amu;
            return m2 * sqrt(erg_lab * (erg_lab + 2 * m1)) /
                   sqrt((m1 + m2) * (m1 + m2) + 2 * m2 * erg_lab) / (hbar * c);
          }()),
          reduced_mass([&erg_cms, &ch, this]() -> real {
            using constants::hbar;
            using constants::c;
            using constants::MeV_per_amu;

            // relativistic-corrected reduced mass [MeV]
            // Eq. 21 of Ingemarsson, 1974
            const auto m1 = ch.proj_mass * MeV_per_amu;
            const auto Ep = m1 + erg_cms;

            return hbar * hbar * c * c * k * k * Ep / (Ep * Ep - m1 * m1);
          }()),
          h2ma([&ch, this]() -> real {
            using constants::hbar;
            using constants::c;

            return hbar * hbar * c * c / (2 * reduced_mass * ch.radius);
          }()),
          sommerfield_param([&ch, this]() -> real {
            using constants::hbar;
            using constants::c;
            using constants::e_sqr;

            // Bohr radius [fm]
            const real ab =
                hbar * hbar * c * c / (reduced_mass * e_sqr * fabs(ch.Zz));

            // dimensionless
            return std::copysign(1., ch.Zz) / (ab * k);
          }()) {}
  };

  /// @brief analytic solutions of free, reduced, radial Schrödinger eqn; also
  /// asymptotic (r → ∞) solutions for non-zero interactions (provided
  /// interaction → 0 faster than 1/r)
  struct Asymptotics {
    /// @brief value of the outgoing asymptotic wavefunction
    cmpl wvfxn_out{};
    /// @brief value of the incoming asymptotic wavefunction
    cmpl wvfxn_in{};
    /// @brief value of the derivative of the outgoing asymptotic wavefunction
    cmpl wvfxn_deriv_out{};
    /// @brief value of the derivative of the incoming asymptotic wavefunction
    cmpl wvfxn_deriv_in{};

    Asymptotics() = default;

    Asymptotics(int l, real k, real r)
        : wvfxn_out(asymptotics::h_plus{l}(k * r)),
          wvfxn_in(asymptotics::h_minus{l}(k * r)),
          wvfxn_deriv_out(asymptotics::d_dz(asymptotics::h_plus{l})(k * r)),
          wvfxn_deriv_in(asymptotics::d_dz(asymptotics::h_minus{l})(k * r)) {}

    static std::vector<Asymptotics> generate_asymptotics(int lmax, real k,
                                                         real r) {
      std::vector<Asymptotics> asym(lmax);
      for (int l = 0; l < lmax; ++l) {
        asym[l] = Asymptotics(l, k, r);
      }
      return asym;
    }
  };

  /// @brief channel energy threshold [MeV]
  real threshold;
  /// @brief channel radius [fm]
  real radius;
  /// @brief target mass [amu]
  real targ_mass;
  /// @brief projectile mass [amu]
  real proj_mass;
  /// @brief product of proton # of target and projectile
  real Zz;
  /// @brief 2s+1, with s projectile spin
  int s;

  Energetics set_erg_cms(real erg_cms) const {
    return Energetics(erg_cms, *this);
  }

  Energetics set_erg_lab(real erg_lab) const {
    using constants::MeV_per_amu;
    const auto m1 = proj_mass * MeV_per_amu;
    const auto m2 = targ_mass * MeV_per_amu;

    const auto erg_cms = erg_lab * m2 / (m1 + m2);
    return set_erg_cms(erg_cms);
  }

  Asymptotics set_angular_momentum(FermionSpinOrbitCoupling am, real k) const {
    return Asymptotics(am.l, k, radius);
  }

  /// @param threshold [Mev]
  /// @param radius [fm]
  /// @param proj_mass [amu]
  /// @param Zp proton number of projectile
  /// @param 2s+1, w/ s spin of projectile
  /// @param targ_mass [amu]
  /// @param Zt proton number of target
  Channel(real threshold, real radius, real proj_mass, int Zp, int s,
          real targ_mass, int Zt)
      : threshold(threshold), radius(radius), targ_mass(targ_mass),
        proj_mass(proj_mass), Zz(Zp * Zt), s(s) {}
};

}; // namespace osiris

#endif
