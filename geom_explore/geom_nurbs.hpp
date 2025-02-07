#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <optional>

namespace geom_examples {

  // build a nurbs curve and evaluate it
  void nurbsExample();

  // build a nurbs curve, do evaluations and add up the time taken
  std::string nurbsPerformanceExample();

}