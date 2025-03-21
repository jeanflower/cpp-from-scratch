#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include "geom.hpp"

namespace geom_examples {

  // Generate N muted colors
  std::vector<uint32_t> generateMutedColors(int count);

  constexpr uint32_t RED     = 0xFF0000;
  constexpr uint32_t GREEN   = 0x00FF00;
  constexpr uint32_t BLUE    = 0x0000FF;
  constexpr uint32_t MAGENTA = 0xFF00FF;
  constexpr uint32_t YELLOW  = 0xFFFF00;
  constexpr uint32_t CYAN    = 0x00FFFF;
  constexpr uint32_t WHITE   = 0xFFFFFF;
  constexpr uint32_t BLACK    = 0x000000;

  const int NUM_SAMPLES = 5000;
  //const int NUM_SAMPLES = 5;

  void clearView();
  void addGeometryToView(const std::vector<PtCollection>& pt_colls);

  void writeGeometryToJSON();

  void addCurveToView(
    const Curve& c,
    double start_t,
    double end_t,
    uint32_t color,
    int displaySize
  );

  void addPointToView(
    const Point& p,
    uint32_t color,
    int displaySize
  );
}