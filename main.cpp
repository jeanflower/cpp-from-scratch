#include <iostream>
#include <string>
#include "geom_explore/geom.hpp"
#include "geom_explore/geomAPI_usage.hpp"
#include "geom_explore/geom_nurbs.hpp"
#include "geom_explore/geom_display.hpp"
#include "boost_explore/boost_usage.hpp"

int main(int argc, char *argv[]) {
  double origin_r = 0, origin_i = 0, radius = 2;
  int num_r = 1000, num_i = 1000;
  char *end = NULL;
  if (argc > 1) {
    origin_r = strtod(argv[1], &end);
    if (argc > 2) {
      origin_i = strtod(argv[2], &end);
      if (argc > 3) {
        radius = strtod(argv[3], &end);
        if (argc > 4) {
          num_i = strtod(argv[4], &end);
          if (argc > 5) {
            num_r = strtod(argv[5], &end);
          }
        }
      }
    }
  }
  std::cout << "Running main() -------------\n";

  //boost_data_types::optionalExample();

  //geomAPI_examples::sphereExample();
  //geomAPI_examples::torusExample();
  //geom_examples::nurbsExample();
  //geom_examples::test_distance_sq_fn();

  //geom_examples::testNewtonRaphson2Dinput();
  //geom_examples::testNewtonRaphson3Dinput();

  // run this with the viewer open to see a plot
  // and interesting for performance comparisons
  geom_examples::fractal(origin_r, origin_i, radius, num_r, num_i);

  // writs any accumulated geometry to a JSON file
  // which the viewer can read
  geom_examples::writeGeometryToJSON();

  std::cout << "Finished running main() -------------\n";
  return 0;
}
