#pragma once

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

  template <typename PtCollection, typename LineCollection>
  void writeGeometryToJSON(
    const PtCollection& pts,  // positions for points
    const LineCollection& polylines // positions for polyline vertices
  );

  void nurbs_example();
}