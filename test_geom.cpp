#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "geom_explore/geom.hpp"
#include "geom_explore/geom_display.hpp"
#include "geom_explore/geom_nurbs.hpp"
#include <iostream>
#include <string>
#include "nr_explore/gaussj.hpp"

TEST_CASE("geom coverage", "[coverage]") {
  geom_examples::nurbsExample();
  const std::string message = geom_examples::nurbsPerformanceExample();
  std::println("{}", message);

  geom_examples::writeGeometryToJSON();

  nr_explore::testGaussj();

  geom_examples::circleExample();

  //REQUIRE(1 + 1 == 3);  // This test will fail
}
TEST_CASE("geom checks", "[checks]") {
  REQUIRE(geom_examples::test_curve_derivs());
}