#include <iostream>
#include <string>
#include "geom_explore/geom.hpp"
//#include "geom_explore/geomAPI_usage.hpp"
//#include "geom_explore/geom_nurbs.hpp"
#include "geom_explore/geom_display.hpp"

int main() {
  std::cout << "Running main() -------------\n";
  //geomAPI_examples::sphereExample();
  //geomAPI_examples::torusExample();
  //geom_examples::nurbsExample();
  geom_examples::test_distance_sq_fn();

  geom_examples::writeGeometryToJSON();

  return 0;
}
