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

  void circleExample();

  bool test_curve_derivs();
}