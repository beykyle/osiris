#include "solver/channel.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

TEST_CASE("Build channel at varying energy") {

  using namespace osiris;
  using namespace osiris::constants;

  constexpr int A = 139;
  constexpr int Z = 54;

  // constants
  constexpr real projectile_mass = n_mass_amu;
  constexpr real target_mass = A;

  constexpr real threshold = 0.; // MeV
  constexpr real ch_radius = 12; // fm

  const auto chn =
      Channel(threshold, ch_radius, projectile_mass, 0, 2, target_mass, Z);
  
  /// p_{3/2}
  const auto am = Channel::FermionSpinOrbitCoupling( 4,  1);

  SECTION("neutron 1MeV") {
    const auto e = chn.set_erg_cms(1.0);
    const auto asym = Channel::Asymptotics(am.l, e.k, chn.radius);
    
    REQUIRE(e.erg_lab == Approx(1.00725658212791));
    REQUIRE(e.k == Approx(0.2189449404));
    REQUIRE(e.h2ma == Approx(1.7374714662639255));
    REQUIRE(e.reduced_mass == Approx(933.7788172169));
    REQUIRE(e.sommerfield_param == Approx(0.));

    REQUIRE(am.l == 1);
    REQUIRE(am.j21 == 4);
    REQUIRE(am.l_dot_s() == 0.5);

    REQUIRE(asym.wvfxn_in.real() == Approx(0.16050015525632866));
    REQUIRE(asym.wvfxn_in.imag() == Approx(-1.0578781142007916));
    
    REQUIRE(asym.wvfxn_out.real() == Approx(0.16050015525632866));
    REQUIRE(asym.wvfxn_out.imag() == Approx(1.0578781142007916));
    
    REQUIRE(asym.wvfxn_deriv_in.real() == Approx(-0.9317486986265555));
    REQUIRE(asym.wvfxn_deriv_in.imag() == Approx(-0.08924255501252532));
    
    REQUIRE(asym.wvfxn_deriv_out.real() == Approx(-0.9317486986265555));
    REQUIRE(asym.wvfxn_deriv_out.imag() == Approx(0.08924255501252532));
  }

  SECTION("neutron 100MeV") {
    const auto e = chn.set_erg_cms(100.0);
    const auto asym = Channel::Asymptotics(am.l, e.k, chn.radius);
    
    REQUIRE(e.erg_lab == Approx(100.725658212791));
    REQUIRE(e.k == Approx(2.24505752889899));
    REQUIRE(e.h2ma == Approx(1.57382520561950));
    REQUIRE(e.reduced_mass == Approx(1030.87383743656));
    REQUIRE(e.sommerfield_param == Approx(0.));
    
    REQUIRE(am.l == 1);
    REQUIRE(am.j21 == 4);
    REQUIRE(am.l_dot_s() == 0.5);
    
    REQUIRE(asym.wvfxn_in.real() == Approx(0.9632903341773514));
    REQUIRE(asym.wvfxn_in.imag() == Approx(-0.2710157201963589));
    
    REQUIRE(asym.wvfxn_out.real() == Approx(0.9632903341773514));
    REQUIRE(asym.wvfxn_out.imag() == Approx(0.2710157201963589));
    
    REQUIRE(asym.wvfxn_deriv_in.real() == Approx(-0.2706920282369296));
    REQUIRE(asym.wvfxn_deriv_in.imag() == Approx(-0.9619511098235556));
    
    REQUIRE(asym.wvfxn_deriv_out.real() == Approx(-0.2706920282369296));
    REQUIRE(asym.wvfxn_deriv_out.imag() == Approx(0.9619511098235556));
  }
}
