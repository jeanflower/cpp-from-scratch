#include <iostream>
#include <string>
#include "geom_explore/geom.hpp"
#include "geom_explore/geomAPI_usage.hpp"
#include "geom_explore/geom_nurbs.hpp"
#include "geom_explore/geom_display.hpp"
#include "boost_usage.hpp"

int main() {
  std::println("Running main() -------------");

  boost_data_types::optionalExample();

  //geomAPI_examples::sphereExample();
  //geomAPI_examples::torusExample();
  //geom_examples::nurbsExample();
  //geom_examples::test_distance_sq_fn();

  //geom_examples::testNewtonRaphson2Dinput();
  //geom_examples::testNewtonRaphson3Dinput();

  // run this with the viewer open to see a plot
  // and interesting for performance comparisons
  //geom_examples::fractal();

  // writs any accumulated geometry to a JSON file
  // which the viewer can read
  //geom_examples::writeGeometryToJSON();

  std::println("Finished running main() -------------");
  return 0;
}
