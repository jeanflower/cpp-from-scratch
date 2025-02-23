#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <tuple>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include "geom.hpp"
#include "geom_display.hpp"
#include <functional>
#include <utility>  // for std::pair


namespace geom_examples {
  // Define a 3D point structure
  Point::Point(double x, double y, double z) : x(x), y(y), z(z) {}

  double Point::X() const { return x; }
  double Point::Y() const { return y; }
  double Point::Z() const { return z; }

  Vector::Vector(double x, double y, double z) : x(x), y(y), z(z) {}

  double Vector::X() const { return x; }
  double Vector::Y() const { return y; }
  double Vector::Z() const { return z; }

  Line::Line(
    Point start,
    Vector direction
  ) : start(start), direction(direction) {}

  Point Line::evaluate(
    double t, 
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{
      start.X() + t * direction.X(), 
      start.Y() + t * direction.Y(), 
      start.Z() + t * direction.Z() 
    };

    if (first_derivative.has_value()) {
      first_derivative->get().x = direction.X();
      first_derivative->get().y = direction.Y();
      first_derivative->get().z = direction.Z();
    }

    if (second_derivative.has_value()) {
      second_derivative->get().x = 0;
      second_derivative->get().y = 0;
      second_derivative->get().z = 0;
    }

    // deliberately break the evaluator to test our tests!!
    //if (t == 1000.0) { // terrible way to pick out a t
    //  p.x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && first_derivative.has_value()) { // terrible way to pick out a t
    //  first_derivative->get().x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && second_derivative.has_value()) { // terrible way to pick out a t
    //  second_derivative->get().x = 88; // nonsense value to return
    //}
    return p;
  }

  Circle::Circle(
    Point origin,
    double r
  ) : origin(origin), radius(r) {}

  Point Circle::evaluate(
    double t, 
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{
      origin.X() + radius * cos(t), 
      origin.Y() + radius * sin(t),
      origin.Z()
    };

    if (first_derivative.has_value()) {
      first_derivative->get().x = -radius * sin(t);
      first_derivative->get().y = radius * cos(t);
    }

    if (second_derivative.has_value()) {
      second_derivative->get().x = -radius * cos(t);
      second_derivative->get().y = -radius * sin(t);
    }

    // deliberately break the evaluator to test our tests!!
    //if (t == 1000.0) { // terrible way to pick out a t
    //  p.x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && first_derivative.has_value()) { // terrible way to pick out a t
    //  first_derivative->get().x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && second_derivative.has_value()) { // terrible way to pick out a t
    //  second_derivative->get().x = 88; // nonsense value to return
    //}
    return p;
  }
  
  void circleExample() {
    Circle c(Point(0.0, 0.0, 0.0), 5.0);
    double t = 1.0;

    Vector first, second;
    Point p = c.evaluate(t, first, second);

    std::cout << "Point: (" << p.x << ", " << p.y << ")\n";
    std::cout << "First Derivative: (" << first.x << ", " << first.y << ")\n";
    std::cout << "Second Derivative: (" << second.x << ", " << second.y << ")\n";
  }

  double calculate_sq_error(
    std::function<std::pair<Vector, Vector>(double)>& geom_evaluator,
    double t,
    const Vector p,
    const Vector& first_deriv,
    double t_step
  ){
    std::pair<Vector, Vector> evals = geom_evaluator(t + t_step); // above t
    Vector p_nearby = evals.first; // reads like position but could actually be a higher deriv
    // std::cout << "p_nearby = " << p_nearby.X() << ", " << p_nearby.Y() << ", " << p_nearby.Z() << "\n";

    Vector estimate( // TODO implement subtraction of Point
      (p_nearby.X() - p.X()) / t_step,
      (p_nearby.Y() - p.Y()) / t_step,
      (p_nearby.Z() - p.Z()) / t_step
    );

    const double sq_error =  // TODO implement subtraction, sq_length of Vector
      (estimate.X() - first_deriv.X()) * (estimate.X() - first_deriv.X()) +
      (estimate.Y() - first_deriv.Y()) * (estimate.Y() - first_deriv.Y()) +
      (estimate.Z() - first_deriv.Z()) * (estimate.Z() - first_deriv.Z());

    std::cout << "estimate with t_step " << t_step << " = "
      << estimate.X() << ", " << estimate.Y() << ", " << estimate.Z()
      << " with sq error " << sq_error
      << "\n";
    return sq_error;
  }

  // can be used for testing derivatives at any order
  // being given evaluators for order n and n+1
  bool test_curve_derivs_generic(
    std::function<std::pair<Vector, Vector>(double)>& geom_evaluator,
    double t
  ) {
    std::pair<Vector, Vector> evals = geom_evaluator(t);
    Vector p = evals.first; // reads like position but could actually be a higher deriv
    Vector first_deriv = evals.second;

    //std::cout << "p = " << p.X() << ", " << p.Y() << ", " << p.Z() << "\n";
    //std::cout << "first_deriv = " << first_deriv.X() << ", " << first_deriv.Y() << ", " << first_deriv.Z() << "\n";

    const double sq_epsilon = 1e-9;

    // convergence means: for all epsilon, there exists a delta such that...
    // for all values of t within delta
    // the error between estimated deriv and reported deriv is less than epsilon
    //
    // we're just picking one epsilon, guessing a delta,
    // and performing a sample of evaluations within delta
    // this isn't quite right in many ways...
    // it assumes a monotonic convergence from our first estimate within tolerance

    // ensure that, if we get close enough, the estimate is close to the reported value
    // to this tolerance
    const double minimal_delta = 1e-8; // if epsilon gets this small, we are close
    const double initial_t_step = 0.001; // start relatively close for a good chance of a good estimate
    const double initial_large_sq_error = std::numeric_limits<double>::max();
    const double t_step_scaling = 0.5; // reduce t_step each time

    // use a parameter-step of t_step to estimate the derivative at t
    double t_step = initial_t_step; 

    double sq_error = initial_large_sq_error;
    bool found_good_estimate = false;
    while (sq_error > sq_epsilon && t_step > minimal_delta) {
      t_step *= t_step_scaling;
      sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      found_good_estimate = sq_error < sq_epsilon;
    }

    if (!found_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, no estimate from above\n";
      return false;
    }

    //std::cout << "found good estimate with delta = " << t_step << "\n";
    bool stayed_good_estimate = true;
    while (stayed_good_estimate && t_step > minimal_delta) {
      t_step *= t_step_scaling;
      const double sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      stayed_good_estimate = sq_error < sq_epsilon;
    }
    if (!stayed_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, lost estimate from above\n";
      return false;
    }
    //std::cout << "kept good estimates to delta = " << t_step << "\n";

    t_step = -initial_t_step; 

    sq_error = initial_large_sq_error;
    found_good_estimate = false;
    while (sq_error > sq_epsilon && -t_step > minimal_delta) {
      t_step *= t_step_scaling;
      sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      found_good_estimate = sq_error < sq_epsilon;
    }

    if (!found_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, no estimate from below\n";
      return false;
    }
    //std::cout << "found good estimate with t_step = " << t_step << "\n";

    stayed_good_estimate = true;
    while (stayed_good_estimate && -t_step > minimal_delta) {
      t_step *= t_step_scaling;
      const double sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      stayed_good_estimate = sq_error < sq_epsilon;
    }
    if (!stayed_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, lost estimate from below\n";
      return false;
    }
    //std::cout << "kept good estimates to delta = " << t_step << "\n";

    return true;
  }

  bool test_curve_derivs_at(
    const Curve& c,
    const double t
  ) {
    bool result = true;
    std::cout << "Test curve derivs at parameter value " << t << "\n";

    std::function<std::pair<Vector, Vector>(double)> position_and_first_deriv_evaluator =
      [&c](double t) -> std::pair<Vector, Vector> {
      Vector firstDeriv;
      Point p = c.evaluate(t, firstDeriv);
      const Vector v(p.X(), p.Y(), p.Z()); // TODO avoid ever copying Point to Vector
      return std::make_pair(
        v, 
        firstDeriv
      );
    };
    result = test_curve_derivs_generic(
      position_and_first_deriv_evaluator,
      t
    );

    if (result) {
      std::function<std::pair<Vector, Vector>(double)> first_and_second_deriv_evaluator =
        [&c](double t) -> std::pair<Vector, Vector> {
        Vector firstDeriv, secondDeriv;
        Point p = c.evaluate(t, firstDeriv, secondDeriv);
        return std::make_pair(
          firstDeriv,
          secondDeriv
        );
      };
      result = test_curve_derivs_generic(
        first_and_second_deriv_evaluator,
        t
      );
    }
    return result;
  }

  bool test_curve_derivs(){

    // select a curve
    const Circle c(Point(0.0, 0.0, 0.0), 12.0);

    // select values of t to analyse for derivative accuracy    
    bool result = true;
    if (result) {
      result = test_curve_derivs_at(c, 0.2);
    }
    if (result) {
      result = test_curve_derivs_at(c, 0);
    }
    if (result) {
      result = test_curve_derivs_at(c, 1000.0);
    }

    return result;
  }

  double DistanceSqFromOrigin::operator() (const double t) {
    Point p = c.evaluate(t);
    return 
      p.X() * p.X() +
      p.Y() * p.Y() +
      p.Z() * p.Z();
  };
  double DistanceSqFromOrigin::df(const double t) {
    Vector deriv;
    Point p = c.evaluate(t, deriv);
    return
      2 * p.X() * deriv.X() +
      2 * p.Y() * deriv.Y() +
      2 * p.Z() * deriv.Z();
  };

  void test_distance_sq_fn() {

    std::function<void(DistanceSqFromOrigin&, double)> print_vals =
      [](DistanceSqFromOrigin& f, double t) -> void {
        const double val = f(t);
        const double dval = f.df(t);
    
        std::cout << "distance from O at " << t << " is " << val << " and df is " << dval << "\n";
      };

    Circle c1(Point(0.0, 0.0, 0.0), 2);
    DistanceSqFromOrigin f1(c1);
  
    print_vals(f1, 0);
    print_vals(f1, 3);

    Circle c2(Point(1.0, 0.0, 0.0), 2);
    DistanceSqFromOrigin f2(c2);
  
    print_vals(f2, 0);
    print_vals(f2, 3);

    Circle c3(Point(1.0, 2.0, 3.0), sqrt(14.0));
    DistanceSqFromOrigin f3(c3);
  
    print_vals(f2, 0);
    print_vals(f2, 3);

    addCurveToView(Line(Point(0,0,0), Vector(1,0,0)), 0.0, 1.0, RED, 1);
    addCurveToView(Line(Point(0,0,0), Vector(0,1,0)), 0.0, 1.0, GREEN, 1);
    addCurveToView(Line(Point(0,0,0), Vector(0,0,1)), 0.0, 1.0, BLUE, 1);

    addCurveToView(c1, 0.0, 2*M_PI, GREEN, 2);
    addCurveToView(c2, 0.0, 2*M_PI, RED, 2);
    addCurveToView(c3, 0.0, 2*M_PI, BLUE, 2);

    // e.g. we might highlight a 'closest point'
    addPointToView(Point(2, 0, 0), YELLOW, 7);

    geom_examples::writeGeometryToJSON();
  }
}
