#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <optional>

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

  struct Vector {
    double x, y, z;

    // Constructor
    Vector(double x = 0, double y = 0, double z = 0);

    // Getter functions
    double X() const;
    double Y() const;
    double Z() const;
  };

  struct PtCollection {
    int color;
    std::vector<Point> pts;
    bool isLine;
  };

  class Curve {
    public:
      // Evaluate curve at parameter t
      // Returns a Point, optionally calculates first and second derivatives
      virtual Point evaluate(
        double t, 
        std::optional<std::reference_wrapper<Vector>> first_derivative = std::nullopt, 
        std::optional<std::reference_wrapper<Vector>> second_derivative = std::nullopt
      ) const = 0;

      virtual ~Curve() = default;
  };

  class Circle : public Curve {
    private:
      double radius;

    public:
      explicit Circle(double r);

      Point evaluate(
        double t, 
        std::optional<std::reference_wrapper<Vector>> first_derivative = std::nullopt, 
        std::optional<std::reference_wrapper<Vector>> second_derivative = std::nullopt
      ) const override;
  };


  constexpr uint32_t RED     = 0xFF0000;
  constexpr uint32_t GREEN   = 0x00FF00;
  constexpr uint32_t BLUE    = 0x0000FF;
  constexpr uint32_t MAGENTA = 0xFF00FF;
  constexpr uint32_t YELLOW  = 0xFFFF00;
  constexpr uint32_t CYAN    = 0x00FFFF;

  const int NUM_SAMPLES = 5000;
  //const int NUM_SAMPLES = 5;

  void addGeometryToView(const std::vector<PtCollection>& pt_colls);

  void writeGeometryToJSON();

  // build a nurbs curve and evaluate it
  void nurbsExample();

  // build a nurbs curve, do evaluations and add up the time taken
  std::string nurbsPerformanceExample();

  void circleExample();

  bool test_curve_derivs();
}