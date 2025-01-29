#pragma once

#include <vector>

namespace geom_examples {

  struct Point {
    double x, y, z;

    // Constructor
    Point(double x = 0, double y = 0, double z = 0);

    // Getter functions
    double X() const;
    double Y() const;
    double Z() const;
  };

  struct PtCollection {
    std::vector<Point> pts;
    int color;
    bool isLine;
  };

  const int NUM_SAMPLES = 5000;
  //const int NUM_SAMPLES = 5;

  void addGeometryToView(const std::vector<PtCollection>& ptColls);

  void writeGeometryToJSON();

  // build a nurbs curve and evaluate it
  void nurbs_example();

  // build a nurbs curve, do evaluations and add up the time taken
  void nurbs_performance_example();
}