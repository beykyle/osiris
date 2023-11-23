#include "extern/legendre_rule.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

TEST_CASE("GL weights shift reduce") {

  SECTION("N=10") {
    const auto quadrature = gl::generate_gauss_legendre_quadrature(10);

    REQUIRE(quadrature[0].weight == Approx(0.03333567));
    REQUIRE(quadrature[3].weight == Approx(0.13463336));
    REQUIRE(quadrature[8].weight == Approx(0.07472567));

    REQUIRE(quadrature[0].abscissa == Approx(0.01304674));
    REQUIRE(quadrature[3].abscissa == Approx(0.2833023));
    REQUIRE(quadrature[8].abscissa == Approx(0.93253168));
  }
}
