#include "potential/potential.hpp"
#include "solver/solver.hpp"
#include "util/constants.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace osiris;

using constants::c;
using constants::hbar;

using Catch::Approx;
/*
using Solver = LagrangeMeshSolver<20>;


TEST_CASE("Yamaguchi analytic s-wave phase shift") {
 // test of R-Matrix solver against analytic potential
 // see Table 16 of 
 // Baye, Daniel. "The Lagrange-mesh method." Physics reports 565 (2015): 1-107.

  // potential
  const auto p = Yamaguchi();

  const auto ch  = Channel(0., 15, constants::n_mass_amu, 0, 2, 
                           constants::p_mass_amu, 1);
  // solver
  Solver solver(ch.radius);

  // S-Wave, 0+
  const auto am   = Channel::AngularMomentum( 2);
  SECTION("0.1MeV") {
    const auto e  = ch.set_erg_cms(0.1);
    const auto k  = e.k;

    const auto swave_k = p.analytic_swave_kmatrix(k);

    REQUIRE(e.erg_cms == 0.1);
    REQUIRE(am.l == 0);
    REQUIRE(k == Approx(sqrt(e.erg_cms/41.472)));
    REQUIRE(e.reduced_mass == Approx(hbar*hbar*c*c/(2*41.472)).epsilon(0.01));
    
    const auto delta = solver.matrices(ch, e, am, p).phase_shift;
    const auto K     = solver.matrices(ch, e, am, p).K();
    
    REQUIRE( delta.real() == Approx(-0.263137 ) );
    REQUIRE( delta.imag() + 1.0 == Approx(1.0) );
    REQUIRE( K.real() == Approx(swave_k) );
    REQUIRE( K.imag() == Approx(0.) );
  }
  
  SECTION("1.89MeV") {
    const auto e  = ch.set_erg_cms(1.89);
    const auto k  = e.k;
    
    const auto swave_k = p.analytic_swave_kmatrix(k);
    
    REQUIRE(am.l == 0);
    REQUIRE( solver.matrices(ch, e, am, p).K().real() == Approx(swave_k) );
  }
  
  SECTION("10MeV") {
    const auto e  = ch.set_erg_cms(10.0);
    const auto k  = e.k;
    
    const auto swave_k = p.analytic_swave_kmatrix(k);

    REQUIRE(e.erg_cms == 10.0);
    REQUIRE(am.l == 0);
    REQUIRE(k == Approx(sqrt(e.erg_cms/41.472)));
    REQUIRE(e.reduced_mass == Approx(hbar*hbar*c*c/(2*41.472)).epsilon(0.01));
    
    const auto delta = solver.matrices(ch, e, am, p).phase_shift;
    const auto K     = solver.matrices(ch, e, am, p).K();
    REQUIRE( delta.real() == Approx(1.4946050256) );
    REQUIRE( delta.imag() + 1.0 == Approx(1.0) );
    REQUIRE( K.real() == Approx(swave_k) );
    REQUIRE( K.imag() == Approx(0.) );
  }
}

*/
