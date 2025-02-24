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
    uint32_t color;
    int displaySize;
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

  class Line : public Curve {
    private:
      Point start;
      Vector direction;

    public:
      explicit Line(
        Point start,
        Vector direction
      );

      Point evaluate(  
        double t, 
        std::optional<std::reference_wrapper<Vector>> first_derivative = std::nullopt, 
        std::optional<std::reference_wrapper<Vector>> second_derivative = std::nullopt
      ) const override;
  };

  class Circle : public Curve {
    private:
      Point center;
      double radius;

    public:
      explicit Circle(Point center, double r);

      Point evaluate(  
        double t, 
        std::optional<std::reference_wrapper<Vector>> first_derivative = std::nullopt, 
        std::optional<std::reference_wrapper<Vector>> second_derivative = std::nullopt
      ) const override;
  };

  void circleExample();

  bool test_curve_derivs();

  class Funcd {
    virtual double operator() (const double x) = 0;
    virtual double df(const double x) = 0;
  };

  class DistanceSqFromOrigin: Funcd {
    private:
      Curve& c;

    public:
      explicit DistanceSqFromOrigin(Circle& c) : c(c) {};

      double operator() (const double x) override;
      double df(const double x) override;
  };

  void test_distance_sq_fn();

  void fractal();
}